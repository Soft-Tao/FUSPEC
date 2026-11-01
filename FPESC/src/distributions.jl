"""Public domain interface for reactant distributions."""
abstract type Distribution end

"""An isotropic Maxwellian distribution with temperature in joules."""
struct Maxwellian <: Distribution
    temperature::Float64

    function Maxwellian(temperature::Real)
        isfinite(temperature) && temperature > 0 || throw(ArgumentError(
            "Maxwellian temperature must be finite and positive",
        ))
        new(Float64(temperature))
    end
end

"""An anisotropic Maxwellian with parallel/perpendicular temperatures."""
struct BiMaxwellian <: Distribution
    parallel_temperature::Float64
    perpendicular_temperature::Float64

    function BiMaxwellian(parallel_temperature::Real,
                          perpendicular_temperature::Real)
        isfinite(parallel_temperature) && parallel_temperature > 0 ||
            throw(ArgumentError(
                "parallel temperature must be finite and positive",
            ))
        isfinite(perpendicular_temperature) && perpendicular_temperature > 0 ||
            throw(ArgumentError(
                "perpendicular temperature must be finite and positive",
            ))
        new(Float64(parallel_temperature), Float64(perpendicular_temperature))
    end
end

"""Distribution tabulated on a strictly increasing (energy, mu) grid."""
struct TabulatedDistribution <: Distribution
    energies::Vector{Float64}
    mus::Vector{Float64}
    values::Matrix{Float64}

    function TabulatedDistribution(energies::AbstractVector{<:Real},
                                   mus::AbstractVector{<:Real},
                                   values::AbstractMatrix{<:Real})
        energy_data = Float64.(energies)
        mu_data = Float64.(mus)
        value_data = Matrix{Float64}(values)
        length(energy_data) >= 2 && all(isfinite, energy_data) &&
            all(diff(energy_data) .> 0) || throw(ArgumentError(
                "tabulated energies must be finite and strictly increasing",
            ))
        length(mu_data) >= 2 && all(isfinite, mu_data) &&
            all(diff(mu_data) .> 0) || throw(ArgumentError(
                "tabulated mus must be finite and strictly increasing",
            ))
        first(mu_data) >= -1.0 && last(mu_data) <= 1.0 ||
            throw(ArgumentError("tabulated mus must lie in [-1, 1]"))
        size(value_data) == (length(energy_data), length(mu_data)) ||
            throw(ArgumentError(
                "tabulated values must have shape (n_energy, n_mu)",
            ))
        all(isfinite, value_data) || throw(ArgumentError(
            "tabulated values must be finite",
        ))
        new(energy_data, mu_data, value_data)
    end
end

@enum DistributionKind::UInt8 begin
    MaxwellianKind = 1
    BiMaxwellianKind = 2
    TabulatedKind = 3
end

"""Type-stable representation consumed by hot kernels."""
struct PackedDistribution
    kind::DistributionKind
    parameters::NTuple{2,Float64}
    energies::Vector{Float64}
    mus::Vector{Float64}
    values::Matrix{Float64}
end

pack_distribution(distribution::Maxwellian) = PackedDistribution(
    MaxwellianKind, (distribution.temperature, 0.0), Float64[], Float64[],
    Matrix{Float64}(undef, 0, 0),
)

pack_distribution(distribution::BiMaxwellian) = PackedDistribution(
    BiMaxwellianKind,
    (distribution.parallel_temperature, distribution.perpendicular_temperature),
    Float64[], Float64[], Matrix{Float64}(undef, 0, 0),
)

pack_distribution(distribution::TabulatedDistribution) = PackedDistribution(
    TabulatedKind, (0.0, 0.0), distribution.energies, distribution.mus,
    distribution.values,
)

@inline function bilinear_value(energies::Vector{Float64},
                                mus::Vector{Float64},
                                values::Matrix{Float64},
                                energy::Float64, mu::Float64)
    if energy < energies[1] || energy > energies[end] ||
       mu < mus[1] || mu > mus[end]
        return 0.0
    end
    energy_index = energy == energies[end] ? length(energies) - 1 :
        searchsortedlast(energies, energy)
    mu_index = mu == mus[end] ? length(mus) - 1 : searchsortedlast(mus, mu)
    energy_fraction = (energy - energies[energy_index]) /
        (energies[energy_index + 1] - energies[energy_index])
    mu_fraction = (mu - mus[mu_index]) /
        (mus[mu_index + 1] - mus[mu_index])
    lower = (1.0 - mu_fraction) * values[energy_index, mu_index] +
        mu_fraction * values[energy_index, mu_index + 1]
    upper = (1.0 - mu_fraction) * values[energy_index + 1, mu_index] +
        mu_fraction * values[energy_index + 1, mu_index + 1]
    return (1.0 - energy_fraction) * lower + energy_fraction * upper
end

@inline function value_vs_momentum(distribution::PackedDistribution,
                                   momentum::Float64, mu::Float64,
                                   mass::Float64)
    if distribution.kind == MaxwellianKind
        temperature = distribution.parameters[1]
        beta = mass / (2.0 * temperature)
        velocity = momentum / mass
        return (beta / pi)^1.5 * velocity^2 *
               exp(-beta * velocity^2) / mass
    elseif distribution.kind == BiMaxwellianKind
        parallel_temperature, perpendicular_temperature =
            distribution.parameters
        beta_parallel = mass / (2.0 * parallel_temperature)
        beta_perpendicular = mass / (2.0 * perpendicular_temperature)
        velocity = momentum / mass
        return sqrt(beta_parallel) * beta_perpendicular / pi^1.5 *
               velocity^2 * exp(
                   -beta_parallel * velocity^2 * mu^2 -
                   beta_perpendicular * velocity^2 * (1.0 - mu^2),
               ) / mass
    end
    energy = 0.5 * momentum^2 / mass
    return bilinear_value(
        distribution.energies, distribution.mus, distribution.values,
        energy, mu,
    ) * momentum / mass
end

"""Uniform importance sample in (|p|, mu, azimuth)."""
@inline function sample_momentum(distribution::PackedDistribution,
                                 mass::Float64,
                                 energy_range::NTuple{2,Float64},
                                 mu_range::NTuple{2,Float64},
                                 seed::UInt64, sample_index::UInt64,
                                 first_dimension::UInt32)
    p_min = sqrt(2.0 * mass * energy_range[1])
    p_max = sqrt(2.0 * mass * energy_range[2])
    magnitude = p_min + counter_uniform(
        seed, sample_index, first_dimension,
    ) * (p_max - p_min)
    mu = mu_range[1] + counter_uniform(
        seed, sample_index, first_dimension + UInt32(1),
    ) * (mu_range[2] - mu_range[1])
    azimuth = 2.0 * pi * counter_uniform(
        seed, sample_index, first_dimension + UInt32(2),
    )
    transverse = sqrt(max(0.0, 1.0 - mu^2))
    momentum = (
        magnitude * transverse * cos(azimuth),
        magnitude * transverse * sin(azimuth),
        magnitude * mu,
    )
    return momentum, value_vs_momentum(distribution, magnitude, mu, mass)
end
