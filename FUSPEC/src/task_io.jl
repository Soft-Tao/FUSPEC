using HDF5

"""A validated prepared HDF5 task and the metadata needed in its result."""
struct PreparedTask{P}
    calculation_kind::Symbol
    problem::P
    R::Vector{Float64}
    Z::Vector{Float64}
    psi_norm::Vector{Float64}
    importance_weight::Vector{Float64}
    B_hat::Union{Nothing,Matrix{Float64}}
    energy_edges::Union{Nothing,Vector{Float64}}
    mu_B::Union{Nothing,Vector{Float64}}
    seed::UInt64
    physics_fingerprint::String
end

function required_object(file, path::AbstractString)
    haskey(file, path) || throw(ArgumentError("missing required HDF5 path /$path"))
    return file[path]
end

read_required(file, path::AbstractString) = read(required_object(file, path))
read_string(file, path::AbstractString) = String(read_required(file, path))

function read_root_attribute(file, name::AbstractString)
    metadata = attributes(file)
    haskey(metadata, name) || throw(ArgumentError(
        "missing required HDF5 root attribute $name",
    ))
    return read(metadata[name])
end

function read_optional_string(file, path::AbstractString)
    haskey(file, path) || return ""
    return String(read(file[path]))
end

function read_float_vector(file, path::AbstractString)
    values = vec(Float64.(read_required(file, path)))
    all(isfinite, values) || throw(ArgumentError("/$path must be finite"))
    return values
end

function validate_zero_offsets(name::AbstractString,
                               offsets::Vector{Int}, count::Int,
                               storage_length::Int)
    length(offsets) == count + 1 || throw(ArgumentError(
        "/$name must have distribution_count + 1 entries",
    ))
    offsets[1] == 0 && offsets[end] == storage_length &&
        all(diff(offsets) .>= 0) || throw(ArgumentError(
            "/$name must be zero-based, non-decreasing, and span storage",
        ))
    return nothing
end

"""Reconstruct a matrix whose flattened HDF5 representation used C order."""
function matrix_from_c_order(values::AbstractVector{<:Real},
                             rows::Int, columns::Int)
    length(values) == rows * columns || throw(DimensionMismatch(
        "C-order flat array length does not match logical matrix shape",
    ))
    return Matrix{Float64}(permutedims(reshape(values, columns, rows)))
end

function flat_c_order(values::AbstractMatrix{<:Real})
    return vec(permutedims(values, (2, 1)))
end

function load_distribution_batch(file, population_path::AbstractString)
    base = "$population_path/distributions"
    kinds = vec(UInt8.(read_required(file, "$base/kind")))
    count = length(kinds)
    parameters_flat = read_float_vector(file, "$base/parameters_flat")
    energy_values = read_float_vector(file, "$base/energy_values")
    mu_values = read_float_vector(file, "$base/mu_values")
    f_values_c = read_float_vector(file, "$base/f_values")
    parameter_offsets = vec(Int.(read_required(file, "$base/parameter_offsets")))
    energy_offsets = vec(Int.(read_required(file, "$base/energy_offsets")))
    mu_offsets = vec(Int.(read_required(file, "$base/mu_offsets")))
    f_offsets = vec(Int.(read_required(file, "$base/f_offsets")))
    validate_zero_offsets(
        "$base/parameter_offsets", parameter_offsets, count,
        length(parameters_flat),
    )
    validate_zero_offsets(
        "$base/energy_offsets", energy_offsets, count, length(energy_values),
    )
    validate_zero_offsets(
        "$base/mu_offsets", mu_offsets, count, length(mu_values),
    )
    validate_zero_offsets(
        "$base/f_offsets", f_offsets, count, length(f_values_c),
    )

    parameters = zeros(Float64, count, 2)
    f_values_julia = Float64[]
    sizehint!(f_values_julia, length(f_values_c))
    for index in 1:count
        kind = kinds[index]
        parameter_start = parameter_offsets[index] + 1
        parameter_stop = parameter_offsets[index + 1]
        parameter_count = parameter_stop - parameter_start + 1
        expected_parameters = kind == UInt8(1) ? 1 : kind == UInt8(2) ? 2 : 0
        parameter_count == expected_parameters || throw(ArgumentError(
            "/$base distribution $(index - 1) has the wrong parameter count",
        ))
        if parameter_count > 0
            parameters[index, 1:parameter_count] .=
                parameters_flat[parameter_start:parameter_stop]
        end

        energy_count = energy_offsets[index + 1] - energy_offsets[index]
        mu_count = mu_offsets[index + 1] - mu_offsets[index]
        f_start = f_offsets[index] + 1
        f_stop = f_offsets[index + 1]
        f_count = f_stop - f_start + 1
        if kind == UInt8(3)
            energy_count >= 2 && mu_count >= 2 &&
                f_count == energy_count * mu_count || throw(ArgumentError(
                    "/$base tabulated distribution $(index - 1) has inconsistent grids",
                ))
            matrix = matrix_from_c_order(
                view(f_values_c, f_start:f_stop), energy_count, mu_count,
            )
            append!(f_values_julia, vec(matrix))
        else
            kind in UInt8.((1, 2)) || throw(ArgumentError(
                "/$base contains unsupported distribution kind $kind",
            ))
            energy_count == 0 && mu_count == 0 && f_count == 0 ||
                throw(ArgumentError(
                    "/$base thermal distribution $(index - 1) owns tabulated data",
                ))
        end
    end
    return PackedDistributionBatch(
        kinds, parameters, energy_offsets .+ 1, mu_offsets .+ 1,
        f_offsets .+ 1, energy_values, mu_values, f_values_julia,
    )
end

function load_population(file, population_path::AbstractString,
                         n_cells::Int)
    density = read_float_vector(file, "$population_path/density")
    length(density) == n_cells || throw(DimensionMismatch(
        "/$population_path/density must have n_cells entries",
    ))
    all(density .>= 0.0) || throw(ArgumentError(
        "/$population_path/density must be non-negative",
    ))
    batch = load_distribution_batch(file, population_path)
    isempty(batch.kinds) && throw(ArgumentError(
        "/$population_path must contain at least one distribution",
    ))
    stored_indices = vec(Int.(read_required(
        file, "$population_path/distribution_index",
    )))
    length(stored_indices) == n_cells || throw(DimensionMismatch(
        "/$population_path/distribution_index must have n_cells entries",
    ))
    indices = Vector{Int}(undef, n_cells)
    for index in eachindex(indices)
        stored = stored_indices[index]
        if stored == -1
            density[index] == 0.0 || throw(ArgumentError(
                "/$population_path cell $(index - 1) has density but no distribution",
            ))
            indices[index] = 1
        elseif 0 <= stored < length(batch.kinds)
            indices[index] = stored + 1
        else
            throw(ArgumentError(
                "/$population_path distribution index $stored is out of range",
            ))
        end
    end
    return density, batch, indices
end

function load_particle(file, role::AbstractString)
    base = "reaction/particles/$role"
    return Particle(
        Symbol(role), Float64(read_required(file, "$base/mass")),
        Int(read_required(file, "$base/charge_number")),
        Int(read_required(file, "$base/mass_number")),
    )
end

function load_reaction(file)
    # The task's differential cross section is always referenced to reactant 1
    # (the projectile), so both products are linked to it. The spectrum kernel
    # mirrors the table with mu -> -mu when the selected product is p4.
    return Reaction(
        :prepared_fusion, load_particle(file, "r1"), load_particle(file, "r2"),
        load_particle(file, "p3"), load_particle(file, "p4");
        link_to_product1=1, link_to_product2=1,
    )
end

function same_particle_definition(left::Particle, right::Particle)
    return left.mass == right.mass &&
        left.charge_number == right.charge_number &&
        left.mass_number == right.mass_number
end

function load_task_reaction(file)
    name = read_string(file, "reaction/name")
    if name == "custom"
        haskey(file, "reaction/particles") || throw(ArgumentError(
            "custom reaction task is missing /reaction/particles",
        ))
        return nothing
    end
    is_builtin_reaction(name) || throw(ArgumentError(
        "unsupported /reaction/name $name",
    ))
    haskey(file, "reaction/particles") && throw(ArgumentError(
        "built-in reaction task must not contain /reaction/particles",
    ))
    haskey(file, "cross_section") && throw(ArgumentError(
        "built-in reaction task must not contain /cross_section",
    ))
    return load_builtin_reaction(name)
end

function load_cells(file)
    R = read_float_vector(file, "cells/R")
    Z = read_float_vector(file, "cells/Z")
    psi_norm = read_float_vector(file, "cells/psi_norm")
    importance = read_float_vector(file, "cells/importance_weight")
    count = length(R)
    length(Z) == count && length(psi_norm) == count &&
        length(importance) == count || throw(DimensionMismatch(
            "all /cells arrays must have the same length",
        ))
    count > 0 || throw(ArgumentError("task contains no cells"))
    all(importance .> 0.0) || throw(ArgumentError(
        "/cells/importance_weight must be positive",
    ))
    return R, Z, psi_norm, importance
end

function load_b_hat(file, n_cells::Int)
    flat = read_float_vector(file, "cells/B_hat_flat")
    matrix = matrix_from_c_order(flat, n_cells, 3)
    for row in axes(matrix, 1)
        isapprox(sum(abs2, view(matrix, row, :)), 1.0; rtol=1e-10, atol=1e-10) ||
            throw(ArgumentError("/cells/B_hat_flat row $(row - 1) is not a unit vector"))
    end
    return matrix
end

function load_ranges(file)
    energy1 = read_float_vector(file, "sampling/reactant1/energy_range")
    energy2 = read_float_vector(file, "sampling/reactant2/energy_range")
    mu1 = read_float_vector(file, "sampling/reactant1/mu_range")
    mu2 = read_float_vector(file, "sampling/reactant2/mu_range")
    all(length(values) == 2 for values in (energy1, energy2, mu1, mu2)) ||
        throw(DimensionMismatch("all sampling ranges must have two entries"))
    return Tuple(energy1), Tuple(energy2), Tuple(mu1), Tuple(mu2)
end

function load_total_cross_section(file)
    read_string(file, "cross_section/kind") == "total-cubic-ppoly" ||
        throw(ArgumentError("reaction-rate task requires total-cubic-ppoly"))
    knots_joule = read_float_vector(file, "cross_section/energy_knots")
    interval_count = length(knots_joule) - 1
    coefficients_si = matrix_from_c_order(
        read_float_vector(file, "cross_section/coefficients_flat"),
        4, interval_count,
    )
    coefficients = copy(coefficients_si)
    for (row, power) in enumerate((3, 2, 1, 0))
        coefficients[row, :] .*= MEV^power / BARN
    end
    return CubicSplineCrossSection(knots_joule ./ MEV, coefficients)
end

function load_differential_cross_section(file)
    read_string(file, "cross_section/kind") ==
        "differential-tensor-bspline" || throw(ArgumentError(
            "spectrum task requires differential-tensor-bspline",
        ))
    degrees = vec(Int.(read_required(file, "cross_section/degrees")))
    degrees == [3, 3] || throw(ArgumentError(
        "only cubic differential tensor B-splines are supported",
    ))
    energy_knots_joule = read_float_vector(file, "cross_section/energy_knots")
    mu_knots = read_float_vector(file, "cross_section/mu_knots")
    rows = length(energy_knots_joule) - 4
    columns = length(mu_knots) - 4
    coefficients = matrix_from_c_order(
        read_float_vector(file, "cross_section/coefficients_flat"), rows, columns,
    ) ./ BARN
    return TensorCubicBSpline(
        energy_knots_joule ./ MEV, mu_knots, coefficients,
    )
end

function normalized_calculation_kind(value::AbstractString)
    value == "reaction-rate-profile" && return :reaction_rate_profile
    value == "product-emission-spectrum" && return :product_emission_spectrum
    throw(ArgumentError("unsupported task calculation_kind $value"))
end

"""Load and validate an FUSPEC prepared task from an HDF5 input bundle."""
function load_task(path::AbstractString)
    isfile(path) || throw(ArgumentError("input task does not exist: $path"))
    return h5open(path, "r") do file
        version = Int(read_root_attribute(file, "schema_version"))
        version == 2 || throw(ArgumentError(
            "unsupported FUSPEC task schema_version $version; expected 2",
        ))
        String(read_root_attribute(file, "unit_system")) == "SI" ||
            throw(ArgumentError("FUSPEC task unit_system must be SI"))
        kind = normalized_calculation_kind(
            read_string(file, "task/calculation_kind"),
        )
        seed = UInt64(read_required(file, "run/seed"))
        fingerprint = read_optional_string(file, "task/physics_fingerprint")
        R, Z, psi_norm, importance = load_cells(file)
        n_cells = length(R)
        resolved = load_task_reaction(file)
        relation = read_string(file, "reaction/population_relation")
        relation in ("same", "distinct") || throw(ArgumentError(
            "/reaction/population_relation must be same or distinct",
        ))
        reaction = resolved === nothing ?
            load_reaction(file) : resolved.reaction
        if resolved !== nothing && resolved.name == "dt-n4he" &&
            relation != "distinct"
            throw(ArgumentError(
                "dt-n4he requires /reaction/population_relation=distinct",
            ))
        end
        if relation == "same" &&
            !same_particle_definition(reaction.reactant1, reaction.reactant2)
            throw(ArgumentError(
                "same-population task requires identical r1 and r2 particles",
            ))
        end
        density1, batch1, indices1 = load_population(
            file, "populations/p0", n_cells,
        )
        if relation == "same"
            density2, batch2, indices2 = density1, batch1, indices1
            pair_factor = 0.5
        else
            density2, batch2, indices2 = load_population(
                file, "populations/p1", n_cells,
            )
            pair_factor = 1.0
        end
        active = (density1 .> 0.0) .& (density2 .> 0.0)
        any(active) || throw(ArgumentError("task has no active reaction cells"))
        cells = build_reaction_rate_cells(
            pair_factor .* density1, density2, importance,
            batch1, batch2, indices1, indices2,
        )
        energy1, energy2, mu1, mu2 = load_ranges(file)

        if kind == :reaction_rate_profile
            haskey(file, "reaction/selected_product") && throw(ArgumentError(
                "reaction-rate task must not select a product",
            ))
            total_cross_section = resolved === nothing ?
                load_total_cross_section(file) : resolved.total_cross_section
            problem = ReactionRateProblem(
                reaction, cells, total_cross_section;
                energy_range1=energy1, energy_range2=energy2,
                mu_range1=mu1, mu_range2=mu2, which_product=1,
            )
            return PreparedTask(
                kind, problem, R, Z, psi_norm, importance, nothing, nothing,
                nothing, seed, fingerprint,
            )
        end

        selected = Int(read_required(file, "reaction/selected_product"))
        selected in (3, 4) || throw(ArgumentError(
            "/reaction/selected_product must be 3 or 4",
        ))
        B_hat = load_b_hat(file, n_cells)
        energy_edges = read_float_vector(file, "spectrum/energy_edges")
        length(energy_edges) >= 2 && all(diff(energy_edges) .> 0.0) ||
            throw(ArgumentError("/spectrum/energy_edges must be increasing"))
        widths = diff(energy_edges)
        all(isapprox(width, widths[1]; rtol=1e-12, atol=0.0) for width in widths) ||
            throw(ArgumentError("spectrum energy bins must be uniform"))
        mu_B = read_float_vector(file, "spectrum/mu_B")
        length(mu_B) > 0 && all(diff(mu_B) .> 0.0) &&
            first(mu_B) >= -1.0 && last(mu_B) <= 1.0 || throw(ArgumentError(
                "/spectrum/mu_B must be increasing inside [-1, 1]",
            ))
        sampling_problem = ReactionRateProblem(
            reaction, cells, ConstantCrossSection(0.0);
            energy_range1=energy1, energy_range2=energy2,
            mu_range1=mu1, mu_range2=mu2,
            which_product=selected == 3 ? 1 : 2,
        )
        differential_cross_section = resolved === nothing ?
            load_differential_cross_section(file) :
            resolved.differential_cross_section
        problem = SpectrumProblem(
            sampling_problem, differential_cross_section;
            emission_mus=mu_B,
            energy_range=(first(energy_edges), last(energy_edges)),
            bins=length(energy_edges) - 1,
        )
        return PreparedTask(
            kind, problem, R, Z, psi_norm, importance, B_hat, energy_edges,
            mu_B, seed, fingerprint,
        )
    end
end

function write_run_metadata(group, task::PreparedTask, config::RunConfig)
    group["calculation_kind"] = replace(String(task.calculation_kind), "_" => "-")
    group["samples_requested"] = config.samples
    group["workers"] = Int64(config.n_workers)
    group["seed"] = config.seed
    group["kinematics"] = "relativistic"
    group["physics_fingerprint"] = task.physics_fingerprint
    return nothing
end

function write_cell_metadata(file, task::PreparedTask)
    cells = create_group(file, "cells")
    cells["R"] = task.R
    cells["Z"] = task.Z
    cells["psi_norm"] = task.psi_norm
    if task.B_hat !== nothing
        cells["B_hat_flat"] = flat_c_order(task.B_hat)
        cells["B_hat_shape"] = Int64[size(task.B_hat)...]
    end
    return nothing
end

function initialize_result_file(file, task::PreparedTask, config::RunConfig)
    metadata = attributes(file)
    metadata["schema_version"] = UInt64(2)
    write_run_metadata(create_group(file, "run"), task, config)
    write_cell_metadata(file, task)
    return nothing
end

"""Write only final normalized output; no accumulator or continuation state."""
function write_result(path::AbstractString,
                      task::PreparedTask{P}, result::EmissionSpectrum,
                      config::RunConfig) where {P<:SpectrumProblem}
    expected = (
        length(task.R), length(task.mu_B), length(task.energy_edges) - 1,
    )
    size(result.accumulator) == expected || throw(DimensionMismatch(
        "spectrum result shape does not match task",
    ))
    h5open(path, "w") do file
        initialize_result_file(file, task, config)
        grids = create_group(file, "grids")
        grids["mu_B"] = task.mu_B
        grids["energy"] = [
            0.5 * (task.energy_edges[index] + task.energy_edges[index + 1])
            for index in 1:(length(task.energy_edges) - 1)
        ]
        result_group = create_group(file, "result")
        n_cells, n_directions, n_bins = expected
        total = n_cells * n_directions * n_bins
        normalized_flat = create_dataset(
            result_group, "normalized_flat", Float64, (total,)
        )
        # HDF5 stores the flattened output in C order (bin, direction, cell).
        # Stream a bounded number of cells at a time to avoid a full-size copy.
        cells_per_chunk = min(n_cells, 64)
        buffer = Vector{Float64}(undef, cells_per_chunk * n_directions * n_bins)
        @inbounds for first_cell in 1:cells_per_chunk:n_cells
            last_cell = min(n_cells, first_cell + cells_per_chunk - 1)
            offset = 1
            for cell_index in first_cell:last_cell
                denominator = result.weight_sums[cell_index]
                for direction_index in 1:n_directions,
                    bin_index in 1:n_bins
                    buffer[offset] = denominator > 0.0 ?
                        result.accumulator[cell_index, direction_index, bin_index] /
                        denominator : 0.0
                    offset += 1
                end
            end
            start_flat = (first_cell - 1) * n_directions * n_bins + 1
            stop_flat = start_flat + offset - 2
            normalized_flat[start_flat:stop_flat] = view(buffer, 1:(offset - 1))
        end
        result_group["normalized_shape"] = Int64[expected...]
    end
    return path
end

function write_result(path::AbstractString,
                      task::PreparedTask{P}, result::ReactionRateResult,
                      config::RunConfig) where {P<:ReactionRateProblem}
    length(result.reaction_rate) == length(task.R) || throw(DimensionMismatch(
        "reaction-rate result length does not match task",
    ))
    physical_rate = [
        result.reaction_rate[index] * task.problem.cells[index].density1 *
            task.problem.cells[index].density2
        for index in eachindex(result.reaction_rate)
    ]
    h5open(path, "w") do file
        initialize_result_file(file, task, config)
        result_group = create_group(file, "result")
        result_group["reaction_rate"] = physical_rate
    end
    return path
end
