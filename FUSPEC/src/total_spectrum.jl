const COUNTER_RNG_ALGORITHM = :fuspec_counter_v1

"""
Inputs for a directional, energy-binned total-spectrum calculation.

`emission_mus` are pitches against the local magnetic field, and the kernel
runs in the field-aligned frame described in section 7 of the
TaskPackageProtocol: the polar axis of every sampled reactant momentum is the
local field, so `emission_directions[index]` is the frame's unit vector whose
polar cosine about that axis is `emission_mus[index]`. The line-of-sight pitch
`n_los . B_hat` is a scalar invariant of the frame and is applied by the
downstream convolution, not here; expressing the emission directions in the
laboratory basis while the sampled momenta stay field-aligned would mix two
reference frames and distort the beam-thermal spectrum.
"""
struct SpectrumProblem{C<:TotalCrossSection,D<:DifferentialCrossSection}
    sampling::ReactionRateProblem{C}
    differential_cross_section::D
    emission_mus::Vector{Float64}
    emission_directions::Vector{Vec3}
    energy_min::Float64
    energy_max::Float64
    bin_width::Float64
    bins::Int
    is_isotropic::Bool
    isotropic_direction_index::Int
end

function SpectrumProblem(sampling::ReactionRateProblem{C},
                         differential_cross_section::D;
                         emission_mus::AbstractVector{<:Real},
                         energy_range::Tuple{<:Real,<:Real},
                         bins::Integer) where
                         {C<:TotalCrossSection,
                          D<:DifferentialCrossSection}
    mus = Float64.(emission_mus)
    isempty(mus) && throw(ArgumentError("emission_mus cannot be empty"))
    all(isfinite, mus) && all(mu -> -1.0 <= mu <= 1.0, mus) ||
        throw(ArgumentError("emission_mus must be finite and lie in [-1, 1]"))
    minimum_energy, maximum_energy = Float64.(energy_range)
    isfinite(minimum_energy) && isfinite(maximum_energy) &&
        minimum_energy < maximum_energy || throw(ArgumentError(
            "energy_range must be finite and increasing",
        ))
    bins > 0 || throw(ArgumentError("bins must be positive"))
    bins <= typemax(Int) || throw(ArgumentError("bins exceed Int"))
    directions = Vec3[
        (sqrt(max(0.0, 1.0 - mu^2)), 0.0, mu) for mu in mus
    ]
    width = (maximum_energy - minimum_energy) / Int(bins)
    isotropic = all(
        cell -> cell.distribution1.kind == MaxwellianKind &&
            cell.distribution2.kind == MaxwellianKind,
        sampling.cells,
    )
    representative_direction = argmax(mus)
    return SpectrumProblem(
        sampling, differential_cross_section, mus, directions,
        minimum_energy, maximum_energy, width, Int(bins), isotropic,
        representative_direction,
    )
end

"""
Authoritative mutable state for deterministic spectrum calculation.

`sample_offset` is the next counter-RNG index for `seed`. `n_samples` counts
all accumulated samples, including samples retained across an explicit seed
change. Normalized values are derived on demand from the raw accumulators and
are not retained in memory.
"""
mutable struct EmissionSpectrum
    accumulator::Array{Float64,3}
    weight_sums::Vector{Float64}
    n_samples::UInt64
    infeasible_samples::UInt64
    seed::Union{Nothing,UInt64}
    sample_offset::UInt64
    rng_algorithm::Symbol
    kernel_version::VersionNumber
end

"""Materialize normalized values for compatibility with the legacy property."""
function normalized_spectrum(result::EmissionSpectrum)
    dimensions = size(result.accumulator)
    normalized = Array{Float64,3}(undef, dimensions)
    Threads.@threads :static for cell_index in axes(result.accumulator, 1)
        @inbounds begin
        denominator = result.weight_sums[cell_index]
        for direction_index in axes(result.accumulator, 2),
            bin_index in axes(result.accumulator, 3)
            normalized[cell_index, direction_index, bin_index] =
                denominator > 0.0 ?
                result.accumulator[cell_index, direction_index, bin_index] /
                denominator : 0.0
        end
        end
    end
    return normalized
end

"""Per-thread scratch for batched direction kinematics."""
struct DirectionBatchScratch
    product_energy::Vector{Float64}
    scattering_cosine::Vector{Float64}
    jacobian::Vector{Float64}
    feasible::Vector{Bool}
end

DirectionBatchScratch(directions::Integer) = DirectionBatchScratch(
    Vector{Float64}(undef, directions),
    Vector{Float64}(undef, directions),
    Vector{Float64}(undef, directions),
    Vector{Bool}(undef, directions),
)

"""Fixed-size per-sample scratch storage for deterministic parallel spectra."""
struct SpectrumWorkspace
    cell_indices::Vector{Int}
    weights::Vector{Float64}
    feasible::Vector{UInt8}
    bin_indices::Matrix{Int}
    contributions::Matrix{Float64}
    direction_scratch::Vector{DirectionBatchScratch}

    function SpectrumWorkspace(problem::SpectrumProblem,
                               block_size::Integer=4096)
        block_size > 0 || throw(ArgumentError(
            "spectrum workspace block_size must be positive",
        ))
        block_size <= typemax(Int) || throw(ArgumentError(
            "spectrum workspace block_size exceeds Int",
        ))
        size = Int(block_size)
        directions = length(problem.emission_directions)
        new(Vector{Int}(undef, size), Vector{Float64}(undef, size),
            Vector{UInt8}(undef, size), Matrix{Int}(undef, directions, size),
            Matrix{Float64}(undef, directions, size),
            [DirectionBatchScratch(directions) for _ in 1:Threads.nthreads()])
    end
end

function EmissionSpectrum(problem::SpectrumProblem)
    dimensions = (
        length(problem.sampling.cells), length(problem.emission_mus),
        problem.bins,
    )
    return EmissionSpectrum(
        zeros(dimensions), zeros(dimensions[1]),
        UInt64(0), UInt64(0), nothing, UInt64(0), COUNTER_RNG_ALGORITHM,
        kernel_version(),
    )
end

function validate_spectrum_result(result::EmissionSpectrum,
                                  problem::SpectrumProblem)
    expected = (
        length(problem.sampling.cells), length(problem.emission_mus),
        problem.bins,
    )
    size(result.accumulator) == expected || throw(DimensionMismatch(
        "result accumulator dimensions do not match the spectrum problem",
    ))
    length(result.weight_sums) == expected[1] || throw(DimensionMismatch(
        "result and problem cell counts differ",
    ))
    result.rng_algorithm == COUNTER_RNG_ALGORITHM || throw(ArgumentError(
        "unsupported RNG algorithm in spectrum state",
    ))
    result.kernel_version == kernel_version() || throw(ArgumentError(
        "spectrum state kernel version does not match this FUSPEC version",
    ))
    return nothing
end

function validate_spectrum_stream(result::EmissionSpectrum,
                                  config::RunConfig)
    result.seed === nothing && return nothing
    if result.seed == config.seed
        result.sample_offset == config.sample_offset || throw(ArgumentError(
            "sample_offset must continue at $(result.sample_offset) for " *
            "seed $(config.seed)",
        ))
    else
        config.sample_offset == 0 || throw(ArgumentError(
            "a new seed must start at sample_offset 0",
        ))
    end
    return nothing
end

"""Replicate the representative direction for isotropic thermal spectra."""
function copy_isotropic_directions!(result::EmissionSpectrum,
                                   problem::SpectrumProblem)
    problem.is_isotropic || return result
    representative = problem.isotropic_direction_index
    @views for direction_index in axes(result.accumulator, 2)
        direction_index == representative && continue
        copyto!(result.accumulator[:, direction_index, :],
                result.accumulator[:, representative, :])
    end
    return result
end

@inline function spectrum_sample_for_cell!(workspace::SpectrumWorkspace,
                                           slot::Int,
                                           problem::SpectrumProblem,
                                           config::RunConfig,
                                           sample_index::UInt64,
                                           cell_index::Int,
                                           worker_index::Int)
    sampling = problem.sampling
    cell = sampling.cells[cell_index]
    momentum1, weight1 = sample_momentum(
        cell.distribution1, sampling.reaction.reactant1.mass,
        sampling.energy_range1, sampling.mu_range1,
        config.seed, sample_index, UInt32(1),
    )
    momentum2, weight2 = sample_momentum(
        cell.distribution2, sampling.reaction.reactant2.mass,
        sampling.energy_range2, sampling.mu_range2,
        config.seed, sample_index, UInt32(4),
    )
    sample_weight = weight1 * weight2
    isfinite(sample_weight) && sample_weight >= 0.0 || throw(ErrorException(
        "total-spectrum sample produced an invalid weight",
    ))
    shared = prepare_kinematics(
        config.kinematics, sampling.reaction, sampling.which_product,
        momentum1, momentum2, config.collective_velocity1,
        config.collective_velocity2,
    )
    workspace.cell_indices[slot] = cell_index
    workspace.weights[slot] = sample_weight
    workspace.feasible[slot] = shared.feasible ? UInt8(1) : UInt8(0)
    energy_cache = prepare_differential_cross_section(
        problem.differential_cross_section, 0.0,
    )
    if shared.feasible
        energy_cache = prepare_differential_cross_section(
            problem.differential_cross_section,
            shared.center_of_mass_energy / MEV,
        )
    end
    scratch = workspace.direction_scratch[worker_index]
    if shared.feasible && !problem.is_isotropic
        direction_kinematics_batch!(
            scratch.product_energy, scratch.scattering_cosine,
            scratch.jacobian, scratch.feasible,
            config.kinematics, shared, problem.emission_directions,
        )
    end
    @inbounds for direction_index in eachindex(problem.emission_directions)
        problem.is_isotropic &&
            direction_index != problem.isotropic_direction_index && continue
        bin_index = 0
        contribution = 0.0
        if shared.feasible
            feasible_flag = false
            energy = 0.0
            cosine = 0.0
            jacobian_value = 0.0
            if problem.is_isotropic
                direction = direction_kinematics(
                    config.kinematics, shared,
                    problem.emission_directions[direction_index],
                )
                feasible_flag = direction.feasible
                energy = direction.product_energy
                cosine = direction.scattering_cosine
                jacobian_value = direction.jacobian
            else
                feasible_flag = scratch.feasible[direction_index]
                energy = scratch.product_energy[direction_index]
                cosine = scratch.scattering_cosine[direction_index]
                jacobian_value = scratch.jacobian[direction_index]
            end
            if feasible_flag
                position = (energy - problem.energy_min) / problem.bin_width
                if 0.0 <= position < problem.bins
                    bin_index = floor(Int, position) + 1
                    mu_query = sampling.which_product == 1 ? cosine : -cosine
                    sigma_barns = evaluate_differential_cross_section(
                        problem.differential_cross_section, energy_cache,
                        mu_query,
                    )
                    isfinite(sigma_barns) || throw(ErrorException(
                        "total-spectrum sample produced a non-finite cross section",
                    ))
                    contribution = cell.density1 * cell.density2 *
                        sigma_barns * BARN * jacobian_value *
                        shared.relative_speed * sample_weight
                end
            end
        end
        workspace.bin_indices[direction_index, slot] = bin_index
        workspace.contributions[direction_index, slot] = contribution
    end
    return nothing
end

@inline function spectrum_sample!(workspace::SpectrumWorkspace, slot::Int,
                                  problem::SpectrumProblem,
                                  config::RunConfig,
                                  sample_index::UInt64)
    cell_index = sample_cell(
        problem.sampling.cell_cdf,
        counter_uniform(config.seed, sample_index, UInt32(0)),
    )
    return spectrum_sample_for_cell!(
        workspace, slot, problem, config, sample_index, cell_index, 1,
    )
end

function fill_spectrum_workspace!(workspace::SpectrumWorkspace,
                                  block_length::Int,
                                  first_sample::UInt64,
                                  problem::SpectrumProblem,
                                  config::RunConfig)
    sampling = problem.sampling
    # Cell selection remains part of the global counter-RNG stream.  After
    # selection, consecutive active-cell ranks are assigned round-robin to
    # logical workers, so a cell and all of its samples have one owner.
    @inbounds for slot in 1:block_length
        sample_index = first_sample + UInt64(slot - 1)
        workspace.cell_indices[slot] = sample_cell(
            sampling.cell_cdf,
            counter_uniform(config.seed, sample_index, UInt32(0)),
        )
    end
    Threads.@threads :static for worker_index in 1:config.n_workers
        @inbounds for slot in 1:block_length
            cell_index = workspace.cell_indices[slot]
            cell_worker(sampling, cell_index, config.n_workers) ==
                worker_index || continue
            spectrum_sample_for_cell!(
                workspace, slot, problem, config,
                first_sample + UInt64(slot - 1), cell_index, worker_index,
            )
        end
    end
    return nothing
end

function validate_spectrum_workspace(workspace::SpectrumWorkspace,
                                     problem::SpectrumProblem)
    block_size = length(workspace.cell_indices)
    block_size > 0 || throw(ArgumentError(
        "spectrum workspace cannot be empty",
    ))
    length(workspace.weights) == block_size &&
        length(workspace.feasible) == block_size &&
        size(workspace.bin_indices) == (
            length(problem.emission_directions), block_size,
        ) && size(workspace.contributions) == size(workspace.bin_indices) ||
        throw(DimensionMismatch(
            "workspace dimensions do not match the spectrum problem",
        ))
    !isempty(workspace.direction_scratch) || throw(DimensionMismatch(
        "workspace direction scratch cannot be empty",
    ))
    all(scratch -> (
        length(scratch.product_energy) == length(problem.emission_directions) &&
        length(scratch.scattering_cosine) == length(problem.emission_directions) &&
        length(scratch.jacobian) == length(problem.emission_directions) &&
        length(scratch.feasible) == length(problem.emission_directions)
    ), workspace.direction_scratch) || throw(DimensionMismatch(
        "workspace direction scratch dimensions do not match the spectrum problem",
    ))
    return block_size
end

function finish_spectrum!(result::EmissionSpectrum, config::RunConfig)
    result.n_samples += config.samples
    result.seed = config.seed
    result.sample_offset = config.sample_offset + config.samples
    return result
end

"""
Run `config.samples` samples and accumulate them into `result`.

Every worker count uses the same block loop below, so the per-sample operation
order does not depend on `julia -t N`. The workspace is sized by the block
length, not by the total sample count.
"""
function calculate_total_spectrum!(result::EmissionSpectrum,
                                   problem::SpectrumProblem,
                                   config::RunConfig)
    block_size = Int(min(max(config.samples, UInt64(1)), UInt64(4096)))
    return calculate_total_spectrum!(
        result, problem, config, SpectrumWorkspace(problem, block_size),
    )
end

function calculate_total_spectrum!(result::EmissionSpectrum,
                                   problem::SpectrumProblem,
                                   config::RunConfig,
                                   workspace::SpectrumWorkspace)
    validate_spectrum_result(result, problem)
    validate_spectrum_stream(result, config)
    config.samples == 0 && return result
    result.n_samples <= typemax(UInt64) - config.samples ||
        throw(ArgumentError("accumulated sample count would overflow UInt64"))
    block_size = validate_spectrum_workspace(workspace, problem)

    next_sample = config.sample_offset
    remaining = config.samples
    while remaining > 0
        block_length = Int(min(remaining, UInt64(block_size)))
        if config.n_workers == 1
            for slot in 1:block_length
                spectrum_sample!(
                    workspace, slot, problem, config,
                    next_sample + UInt64(slot - 1),
                )
            end
        else
            fill_spectrum_workspace!(
                workspace, block_length, next_sample, problem, config,
            )
        end
        @inbounds for slot in 1:block_length
            cell_index = workspace.cell_indices[slot]
            result.weight_sums[cell_index] += workspace.weights[slot]
            if workspace.feasible[slot] == UInt8(1)
                @inbounds for direction_index in eachindex(problem.emission_directions)
                    problem.is_isotropic &&
                        direction_index != problem.isotropic_direction_index && continue
                    bin_index = workspace.bin_indices[direction_index, slot]
                    if bin_index > 0
                        result.accumulator[
                            cell_index, direction_index, bin_index,
                        ] += workspace.contributions[direction_index, slot]
                    end
                end
            else
                result.infeasible_samples += UInt64(1)
            end
        end
        next_sample += UInt64(block_length)
        remaining -= UInt64(block_length)
    end
    finish_spectrum!(result, config)
    return copy_isotropic_directions!(result, problem)
end

function calculate_total_spectrum(problem::SpectrumProblem,
                                  config::RunConfig)
    result = EmissionSpectrum(problem)
    return calculate_total_spectrum!(result, problem, config)
end
