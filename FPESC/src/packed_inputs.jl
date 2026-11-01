"""
Columnar distribution input used by language adapters and bulk I/O.

Kinds are 1=Maxwellian, 2=BiMaxwellian and 3=tabulated. Offsets are Julia
one-based, half-open indices into the flattened tabulated arrays. Values for
each tabulated distribution are stored in Julia column-major order.
"""
struct PackedDistributionBatch
    kinds::Vector{UInt8}
    parameters::Matrix{Float64}
    energy_offsets::Vector{Int}
    mu_offsets::Vector{Int}
    value_offsets::Vector{Int}
    energies::Vector{Float64}
    mus::Vector{Float64}
    values::Vector{Float64}

    function PackedDistributionBatch(
        kinds::Vector{UInt8}, parameters::Matrix{Float64},
        energy_offsets::Vector{Int}, mu_offsets::Vector{Int},
        value_offsets::Vector{Int}, energies::Vector{Float64},
        mus::Vector{Float64}, values::Vector{Float64}, ::Val{:validated},
    )
        new(kinds, parameters, energy_offsets, mu_offsets, value_offsets,
            energies, mus, values)
    end
end

function validate_batch_offsets(name::AbstractString, offsets::Vector{Int},
                                count::Int, storage_length::Int)
    length(offsets) == count + 1 || throw(ArgumentError(
        "$name offsets must have length distribution_count + 1",
    ))
    offsets[1] == 1 && offsets[end] == storage_length + 1 &&
        all(diff(offsets) .>= 0) || throw(ArgumentError(
            "$name offsets must be one-based, non-decreasing and span storage",
        ))
    return nothing
end

function PackedDistributionBatch(kinds::AbstractVector{<:Integer},
                                 parameters::AbstractMatrix{<:Real},
                                 energy_offsets::AbstractVector{<:Integer},
                                 mu_offsets::AbstractVector{<:Integer},
                                 value_offsets::AbstractVector{<:Integer},
                                 energies::AbstractVector{<:Real},
                                 mus::AbstractVector{<:Real},
                                 values::AbstractVector{<:Real})
    kind_data = UInt8.(kinds)
    parameter_data = Matrix{Float64}(parameters)
    energy_offset_data = Int.(energy_offsets)
    mu_offset_data = Int.(mu_offsets)
    value_offset_data = Int.(value_offsets)
    energy_data = Float64.(energies)
    mu_data = Float64.(mus)
    value_data = Float64.(values)
    count = length(kind_data)
    size(parameter_data) == (count, 2) || throw(ArgumentError(
        "distribution parameters must have shape (count, 2)",
    ))
    all(isfinite, parameter_data) || throw(ArgumentError(
        "distribution parameters must be finite",
    ))
    validate_batch_offsets(
        "energy", energy_offset_data, count, length(energy_data),
    )
    validate_batch_offsets("mu", mu_offset_data, count, length(mu_data))
    validate_batch_offsets(
        "value", value_offset_data, count, length(value_data),
    )
    all(isfinite, energy_data) && all(isfinite, mu_data) &&
        all(isfinite, value_data) || throw(ArgumentError(
            "tabulated distribution storage must be finite",
        ))
    for index in 1:count
        kind = kind_data[index]
        kind in UInt8.((1, 2, 3)) || throw(ArgumentError(
            "unknown distribution kind $kind at index $index",
        ))
        energy_count = energy_offset_data[index + 1] -
            energy_offset_data[index]
        mu_count = mu_offset_data[index + 1] - mu_offset_data[index]
        value_count = value_offset_data[index + 1] -
            value_offset_data[index]
        if kind == UInt8(1)
            parameter_data[index, 1] > 0.0 || throw(ArgumentError(
                "Maxwellian temperature must be positive",
            ))
            energy_count == mu_count == value_count == 0 ||
                throw(ArgumentError("thermal distributions cannot own grids"))
        elseif kind == UInt8(2)
            parameter_data[index, 1] > 0.0 &&
                parameter_data[index, 2] > 0.0 || throw(ArgumentError(
                    "BiMaxwellian temperatures must be positive",
                ))
            energy_count == mu_count == value_count == 0 ||
                throw(ArgumentError("thermal distributions cannot own grids"))
        else
            energy_count >= 2 && mu_count >= 2 &&
                value_count == energy_count * mu_count || throw(ArgumentError(
                    "tabulated distribution grid and value spans disagree",
                ))
        end
    end
    return PackedDistributionBatch(
        kind_data, parameter_data, energy_offset_data, mu_offset_data,
        value_offset_data, energy_data, mu_data, value_data, Val(:validated),
    )
end

function distribution_at(batch::PackedDistributionBatch, index::Int)
    kind = batch.kinds[index]
    kind == UInt8(1) && return Maxwellian(batch.parameters[index, 1])
    kind == UInt8(2) && return BiMaxwellian(
        batch.parameters[index, 1], batch.parameters[index, 2],
    )
    energy_range = batch.energy_offsets[index]:(
        batch.energy_offsets[index + 1] - 1
    )
    mu_range = batch.mu_offsets[index]:(batch.mu_offsets[index + 1] - 1)
    value_range = batch.value_offsets[index]:(
        batch.value_offsets[index + 1] - 1
    )
    energy_count = length(energy_range)
    mu_count = length(mu_range)
    return TabulatedDistribution(
        view(batch.energies, energy_range), view(batch.mus, mu_range),
        reshape(view(batch.values, value_range), energy_count, mu_count),
    )
end

"""Build all reaction-rate cells in one adapter boundary crossing."""
function build_reaction_rate_cells(density1::AbstractVector{<:Real},
                                   density2::AbstractVector{<:Real},
                                   weights::AbstractVector{<:Real},
                                   batch1::PackedDistributionBatch,
                                   batch2::PackedDistributionBatch)
    count = length(density1)
    length(density2) == count && length(weights) == count &&
        length(batch1.kinds) == count && length(batch2.kinds) == count ||
        throw(DimensionMismatch("packed cell inputs have different lengths"))
    return build_reaction_rate_cells(
        density1, density2, weights, batch1, batch2, collect(1:count),
        collect(1:count),
    )
end

"""Build cells while sharing repeated packed distribution storage."""
function build_reaction_rate_cells(density1::AbstractVector{<:Real},
                                   density2::AbstractVector{<:Real},
                                   weights::AbstractVector{<:Real},
                                   batch1::PackedDistributionBatch,
                                   batch2::PackedDistributionBatch,
                                   distribution_indices1::AbstractVector{<:Integer},
                                   distribution_indices2::AbstractVector{<:Integer})
    count = length(density1)
    length(density2) == count && length(weights) == count &&
        length(distribution_indices1) == count &&
        length(distribution_indices2) == count || throw(DimensionMismatch(
            "packed cell inputs have different lengths",
        ))
    all(index -> 1 <= index <= length(batch1.kinds),
        distribution_indices1) || throw(ArgumentError(
            "reactant-1 distribution index is out of bounds",
        ))
    all(index -> 1 <= index <= length(batch2.kinds),
        distribution_indices2) || throw(ArgumentError(
            "reactant-2 distribution index is out of bounds",
        ))
    distributions1 = [
        pack_distribution(distribution_at(batch1, index))
        for index in eachindex(batch1.kinds)
    ]
    distributions2 = [
        pack_distribution(distribution_at(batch2, index))
        for index in eachindex(batch2.kinds)
    ]
    cells = Vector{ReactionRateCell}(undef, count)
    for index in 1:count
        cells[index] = ReactionRateCell(
            Float64(density1[index]), Float64(density2[index]),
            Float64(weights[index]),
            distributions1[distribution_indices1[index]],
            distributions2[distribution_indices2[index]],
        )
    end
    return cells
end
