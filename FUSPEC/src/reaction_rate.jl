const BARN = 1.0e-28

"""A reaction-rate cell with type-stable packed reactant distributions."""
struct ReactionRateCell
    density1::Float64
    density2::Float64
    weight::Float64
    distribution1::PackedDistribution
    distribution2::PackedDistribution

    function ReactionRateCell(density1::Float64, density2::Float64,
                              weight::Float64,
                              distribution1::PackedDistribution,
                              distribution2::PackedDistribution)
        isfinite(density1) && density1 >= 0 || throw(ArgumentError(
            "density1 must be finite and non-negative",
        ))
        isfinite(density2) && density2 >= 0 || throw(ArgumentError(
            "density2 must be finite and non-negative",
        ))
        isfinite(weight) && weight >= 0 || throw(ArgumentError(
            "cell weight must be finite and non-negative",
        ))
        new(density1, density2, weight, distribution1, distribution2)
    end
end

"""Packed inputs for the reaction-rate kernel."""
struct ReactionRateProblem{C<:TotalCrossSection}
    reaction::Reaction
    cells::Vector{ReactionRateCell}
    cross_section::C
    energy_range1::NTuple{2,Float64}
    energy_range2::NTuple{2,Float64}
    mu_range1::NTuple{2,Float64}
    mu_range2::NTuple{2,Float64}
    cell_cdf::Vector{Float64}
    active_cell_ranks::Vector{Int}
    which_product::Int8
end

function validate_range(name::AbstractString, values::Tuple{<:Real,<:Real};
                        lower_bound::Float64=0.0,
                        upper_bound::Float64=Inf)
    result = (Float64(values[1]), Float64(values[2]))
    all(isfinite, result) || throw(ArgumentError("$name must be finite"))
    lower_bound <= result[1] < result[2] <= upper_bound ||
        throw(ArgumentError("$name is outside its valid range"))
    return result
end

function build_cell_cdf(cells::Vector{ReactionRateCell})
    isempty(cells) && throw(ArgumentError("reaction-rate cells cannot be empty"))
    cdf = zeros(Float64, length(cells) + 1)
    for (index, cell) in pairs(cells)
        contribution = cell.weight * cell.density1 * cell.density2
        isfinite(contribution) || throw(ArgumentError(
            "cell sampling contribution must be finite",
        ))
        cdf[index + 1] = cdf[index] + contribution
    end
    cdf[end] > 0 || throw(ArgumentError(
        "reaction-rate problem has no active cells",
    ))
    cdf ./= cdf[end]
    cdf[end] = 1.0
    return cdf
end

"""Assign consecutive ranks to active cells; inactive cells retain rank zero."""
function build_active_cell_ranks(cdf::Vector{Float64})
    ranks = zeros(Int, length(cdf) - 1)
    active_rank = 0
    for cell_index in eachindex(ranks)
        if cdf[cell_index + 1] > cdf[cell_index]
            active_rank += 1
            ranks[cell_index] = active_rank
        end
    end
    return ranks
end

"""Logical worker that exclusively owns an active cell during sampling."""
@inline function cell_worker(problem::ReactionRateProblem,
                             cell_index::Int, n_workers::Int)
    active_rank = problem.active_cell_ranks[cell_index]
    return active_rank == 0 ? 0 : mod1(active_rank, n_workers)
end

function ReactionRateProblem(reaction::Reaction,
                             cells::AbstractVector{ReactionRateCell},
                             cross_section::C;
                             energy_range1::Tuple{<:Real,<:Real},
                             energy_range2::Tuple{<:Real,<:Real},
                             mu_range1::Tuple{<:Real,<:Real}=(-1.0, 1.0),
                             mu_range2::Tuple{<:Real,<:Real}=(-1.0, 1.0),
                             which_product::Integer=1) where
                             {C<:TotalCrossSection}
    owned_cells = collect(cells)
    range1 = validate_range("energy_range1", energy_range1)
    range2 = validate_range("energy_range2", energy_range2)
    pitch1 = validate_range(
        "mu_range1", mu_range1; lower_bound=-1.0, upper_bound=1.0,
    )
    pitch2 = validate_range(
        "mu_range2", mu_range2; lower_bound=-1.0, upper_bound=1.0,
    )
    which_product in (1, 2) || throw(ArgumentError(
        "which_product must be 1 or 2",
    ))
    cell_cdf = build_cell_cdf(owned_cells)
    return ReactionRateProblem(
        reaction, owned_cells, cross_section, range1, range2, pitch1, pitch2,
        cell_cdf, build_active_cell_ranks(cell_cdf), Int8(which_product),
    )
end

"""Mutable accumulators and counters for reaction-rate calculation."""
mutable struct ReactionRateResult
    rate_accum::Vector{Float64}
    weight_sums::Vector{Float64}
    reaction_rate::Vector{Float64}
    samples_completed::UInt64
    infeasible_samples::UInt64
end

ReactionRateResult(n_cells::Integer) = ReactionRateResult(
    zeros(n_cells), zeros(n_cells), zeros(n_cells), UInt64(0), UInt64(0),
)

"""Fixed-size scratch storage for deterministic parallel reaction rates."""
struct ReactionRateWorkspace
    cell_indices::Vector{Int}
    weights::Vector{Float64}
    rates::Vector{Float64}
    feasible::Vector{UInt8}

    function ReactionRateWorkspace(block_size::Integer=4096)
        block_size > 0 || throw(ArgumentError(
            "reaction-rate workspace block_size must be positive",
        ))
        block_size <= typemax(Int) || throw(ArgumentError(
            "reaction-rate workspace block_size exceeds Int",
        ))
        size = Int(block_size)
        new(Vector{Int}(undef, size), Vector{Float64}(undef, size),
            Vector{Float64}(undef, size), Vector{UInt8}(undef, size))
    end
end

@inline function sample_cell(cdf::Vector{Float64}, uniform::Float64)
    return min(searchsortedlast(cdf, uniform), length(cdf) - 1)
end

function validate_reaction_rate_result(result::ReactionRateResult,
                                       problem::ReactionRateProblem,
                                       config::RunConfig)
    n_cells = length(problem.cells)
    length(result.rate_accum) == n_cells &&
        length(result.weight_sums) == n_cells &&
        length(result.reaction_rate) == n_cells || throw(DimensionMismatch(
            "result and problem cell counts differ",
        ))
    result.samples_completed <= typemax(UInt64) - config.samples ||
        throw(ArgumentError("accumulated sample count would overflow UInt64"))
    return nothing
end

@inline function reaction_rate_sample_for_cell(problem::ReactionRateProblem,
                                               config::RunConfig,
                                               sample_index::UInt64,
                                               cell_index::Int)
    cell = problem.cells[cell_index]
    momentum1, weight1 = sample_momentum(
        cell.distribution1, problem.reaction.reactant1.mass,
        problem.energy_range1, problem.mu_range1,
        config.seed, sample_index, UInt32(1),
    )
    momentum2, weight2 = sample_momentum(
        cell.distribution2, problem.reaction.reactant2.mass,
        problem.energy_range2, problem.mu_range2,
        config.seed, sample_index, UInt32(4),
    )
    sample_weight = weight1 * weight2
    isfinite(sample_weight) && sample_weight >= 0 || throw(ErrorException(
        "reaction-rate sample produced an invalid weight",
    ))
    shared = prepare_kinematics(
        config.kinematics, problem.reaction, problem.which_product,
        momentum1, momentum2, config.collective_velocity1,
        config.collective_velocity2,
    )
    if shared.feasible
        sigma_barns = cross_section_barns(
            problem.cross_section, shared.center_of_mass_energy / MEV,
        )
        isfinite(sigma_barns) || throw(ErrorException(
            "reaction-rate sample produced a non-finite cross section",
        ))
        rate = sigma_barns * BARN * shared.relative_speed * sample_weight
        return sample_weight, rate, UInt8(1)
    end
    return sample_weight, 0.0, UInt8(0)
end

@inline function reaction_rate_sample(problem::ReactionRateProblem,
                                      config::RunConfig,
                                      sample_index::UInt64)
    cell_index = sample_cell(
        problem.cell_cdf,
        counter_uniform(config.seed, sample_index, UInt32(0)),
    )
    weight, rate, feasible = reaction_rate_sample_for_cell(
        problem, config, sample_index, cell_index,
    )
    return cell_index, weight, rate, feasible
end

function fill_reaction_rate_workspace!(workspace::ReactionRateWorkspace,
                                       block_length::Int,
                                       first_sample::UInt64,
                                       problem::ReactionRateProblem,
                                       config::RunConfig)
    # Select cells first, then give every active cell to exactly one logical
    # worker.  Active-cell ranks are consecutive, so round-robin ownership
    # differs by at most one cell across workers.  A worker handles every
    # sample belonging to its cells; the caller subsequently reduces the
    # per-sample workspace in counter-RNG sample order.
    for slot in 1:block_length
        sample_index = first_sample + UInt64(slot - 1)
        workspace.cell_indices[slot] = sample_cell(
            problem.cell_cdf,
            counter_uniform(config.seed, sample_index, UInt32(0)),
        )
    end
    Threads.@threads :static for worker_index in 1:config.n_workers
        for slot in 1:block_length
            cell_index = workspace.cell_indices[slot]
            cell_worker(problem, cell_index, config.n_workers) ==
                worker_index || continue
            sample_index = first_sample + UInt64(slot - 1)
            weight, rate, feasible = reaction_rate_sample_for_cell(
                problem, config, sample_index, cell_index,
            )
            workspace.weights[slot] = weight
            workspace.rates[slot] = rate
            workspace.feasible[slot] = feasible
        end
    end
    return nothing
end

function normalize_reaction_rate!(result::ReactionRateResult)
    for cell_index in eachindex(result.reaction_rate)
        denominator = result.weight_sums[cell_index]
        result.reaction_rate[cell_index] = denominator > 0 ?
            result.rate_accum[cell_index] / denominator : 0.0
    end
    return result
end

"""
Run `config.samples` samples and accumulate them into `result`.

Every worker count uses the same block loop below, so the reduction order does
not depend on `julia -t N`.
"""
function calculate_reaction_rate!(result::ReactionRateResult,
                                  problem::ReactionRateProblem,
                                  config::RunConfig)
    block_size = Int(min(max(config.samples, UInt64(1)), UInt64(4096)))
    return calculate_reaction_rate!(
        result, problem, config, ReactionRateWorkspace(block_size),
    )
end

function calculate_reaction_rate!(result::ReactionRateResult,
                                  problem::ReactionRateProblem,
                                  config::RunConfig,
                                  workspace::ReactionRateWorkspace)
    validate_reaction_rate_result(result, problem, config)
    config.samples == 0 && return result
    workspace_size = length(workspace.cell_indices)
    workspace_size > 0 || throw(ArgumentError(
        "reaction-rate workspace cannot be empty",
    ))

    next_sample = config.sample_offset
    remaining = config.samples
    while remaining > 0
        block_length = Int(min(remaining, UInt64(workspace_size)))
        if config.n_workers == 1
            for slot in 1:block_length
                sample_index = next_sample + UInt64(slot - 1)
                cell_index, weight, rate, feasible = reaction_rate_sample(
                    problem, config, sample_index,
                )
                workspace.cell_indices[slot] = cell_index
                workspace.weights[slot] = weight
                workspace.rates[slot] = rate
                workspace.feasible[slot] = feasible
            end
        else
            fill_reaction_rate_workspace!(
                workspace, block_length, next_sample, problem, config,
            )
        end
        for slot in 1:block_length
            cell_index = workspace.cell_indices[slot]
            result.weight_sums[cell_index] += workspace.weights[slot]
            if workspace.feasible[slot] == UInt8(1)
                result.rate_accum[cell_index] += workspace.rates[slot]
            else
                result.infeasible_samples += UInt64(1)
            end
            result.samples_completed += UInt64(1)
        end
        next_sample += UInt64(block_length)
        remaining -= UInt64(block_length)
    end
    return normalize_reaction_rate!(result)
end

function calculate_reaction_rate(problem::ReactionRateProblem,
                                 config::RunConfig)
    result = ReactionRateResult(length(problem.cells))
    return calculate_reaction_rate!(result, problem, config)
end
