"""Interface for total cross sections evaluated in barns."""
abstract type TotalCrossSection end

"""Interface for direction-dependent cross sections evaluated in barns/sr."""
abstract type DifferentialCrossSection end

"""Constant total cross section, useful for formula tests."""
struct ConstantCrossSection <: TotalCrossSection
    barns::Float64

    function ConstantCrossSection(barns::Real)
        isfinite(barns) && barns >= 0 || throw(ArgumentError(
            "cross section must be finite and non-negative",
        ))
        new(Float64(barns))
    end
end

@inline cross_section_barns(model::ConstantCrossSection, energy_mev::Float64) =
    model.barns

"""One energy interval of the Bosch-Hale parameterization."""
struct BoschHaleBranch
    energy_min_kev::Float64
    energy_max_kev::Float64
    bg_sqrt_kev::Float64
    A::NTuple{5,Float64}
    B::NTuple{4,Float64}

    function BoschHaleBranch(energy_min_kev::Real, energy_max_kev::Real,
                             bg_sqrt_kev::Real,
                             A::AbstractVector{<:Real},
                             B::AbstractVector{<:Real})
        minimum = Float64(energy_min_kev)
        maximum = Float64(energy_max_kev)
        bg = Float64(bg_sqrt_kev)
        isfinite(minimum) && isfinite(maximum) && 0.0 <= minimum < maximum ||
            throw(ArgumentError(
                "Bosch-Hale energy range must be finite and increasing",
            ))
        isfinite(bg) && bg > 0.0 || throw(ArgumentError(
            "Bosch-Hale BG must be finite and positive",
        ))
        length(A) == 5 && all(isfinite, A) || throw(ArgumentError(
            "Bosch-Hale A must contain five finite coefficients",
        ))
        length(B) == 4 && all(isfinite, B) || throw(ArgumentError(
            "Bosch-Hale B must contain four finite coefficients",
        ))
        new(
            minimum, maximum, bg,
            (Float64(A[1]), Float64(A[2]), Float64(A[3]),
             Float64(A[4]), Float64(A[5])),
            (Float64(B[1]), Float64(B[2]), Float64(B[3]), Float64(B[4])),
        )
    end
end

"""Bosch-Hale total cross section with one or more energy branches."""
struct BoschHaleCrossSection <: TotalCrossSection
    branches::Vector{BoschHaleBranch}

    function BoschHaleCrossSection(branches::AbstractVector{BoschHaleBranch})
        isempty(branches) && throw(ArgumentError(
            "Bosch-Hale cross section requires at least one branch",
        ))
        owned = collect(branches)
        for index in 2:length(owned)
            owned[index].energy_min_kev >= owned[index - 1].energy_max_kev ||
                throw(ArgumentError(
                    "Bosch-Hale branches must be ordered without overlap",
                ))
        end
        new(owned)
    end
end

@inline function horner(coefficients::NTuple{N,Float64},
                        energy::Float64) where {N}
    result = 0.0
    @inbounds for index in N:-1:1
        result = coefficients[index] + energy * result
    end
    return result
end

@inline function cross_section_barns(model::BoschHaleCrossSection,
                                     energy_mev::Float64)
    energy_mev > 0.0 || return 0.0
    energy_kev = energy_mev * 1000.0
    branch = model.branches[end]
    @inbounds for candidate in model.branches
        if energy_kev <= candidate.energy_max_kev
            branch = candidate
            break
        end
    end
    s_factor = horner(branch.A, energy_kev) /
        (1.0 + energy_kev * horner(branch.B, energy_kev))
    sigma_mb = s_factor /
        (energy_kev * exp(branch.bg_sqrt_kev / sqrt(energy_kev)))
    return isfinite(sigma_mb) ? max(0.0, sigma_mb) / 1000.0 : sigma_mb
end

"""
Piecewise cubic power-basis spline compatible with SciPy CubicSpline output.

For interval `i`, coefficients are ordered cubic, quadratic, linear, constant
and evaluated at `dx = energy_mev - knots[i]`. Values outside the knot range
use the first or final polynomial, matching the stage-0 extrapolation contract.
"""
struct CubicSplineCrossSection <: TotalCrossSection
    knots::Vector{Float64}
    coefficients::Matrix{Float64}

    function CubicSplineCrossSection(knots::AbstractVector{<:Real},
                                     coefficients::AbstractMatrix{<:Real})
        knots64 = Float64.(knots)
        coefficients64 = Matrix{Float64}(coefficients)
        length(knots64) >= 2 || throw(ArgumentError(
            "a cubic spline needs at least two knots",
        ))
        all(isfinite, knots64) && all(diff(knots64) .> 0) ||
            throw(ArgumentError("spline knots must be finite and increasing"))
        size(coefficients64) == (4, length(knots64) - 1) ||
            throw(ArgumentError(
                "spline coefficients must have shape (4, length(knots)-1)",
            ))
        all(isfinite, coefficients64) || throw(ArgumentError(
            "spline coefficients must be finite",
        ))
        new(knots64, coefficients64)
    end
end

@inline function cross_section_barns(model::CubicSplineCrossSection,
                                     energy_mev::Float64)
    interval = if energy_mev < model.knots[1]
        1
    elseif energy_mev >= model.knots[end]
        length(model.knots) - 1
    else
        searchsortedlast(model.knots, energy_mev)
    end
    dx = energy_mev - model.knots[interval]
    cubic = model.coefficients[1, interval]
    quadratic = model.coefficients[2, interval]
    linear = model.coefficients[3, interval]
    constant = model.coefficients[4, interval]
    return ((cubic * dx + quadratic) * dx + linear) * dx + constant
end

"""Tensor-product cubic B-spline used by custom differential cross sections."""
struct TensorCubicBSpline <: DifferentialCrossSection
    energy_knots::Vector{Float64}
    mu_knots::Vector{Float64}
    coefficients::Matrix{Float64}

    function TensorCubicBSpline(energy_knots::AbstractVector{<:Real},
                                mu_knots::AbstractVector{<:Real},
                                coefficients::AbstractMatrix{<:Real})
        energy_data = Float64.(energy_knots)
        mu_data = Float64.(mu_knots)
        coefficient_data = Matrix{Float64}(coefficients)
        length(energy_data) >= 8 && all(isfinite, energy_data) &&
            all(diff(energy_data) .>= 0.0) || throw(ArgumentError(
                "energy knots must be finite and non-decreasing",
            ))
        length(mu_data) >= 8 && all(isfinite, mu_data) &&
            all(diff(mu_data) .>= 0.0) || throw(ArgumentError(
                "mu knots must be finite and non-decreasing",
            ))
        expected_shape = (length(energy_data) - 4, length(mu_data) - 4)
        size(coefficient_data) == expected_shape || throw(ArgumentError(
            "B-spline coefficient shape does not match knot counts",
        ))
        all(isfinite, coefficient_data) || throw(ArgumentError(
            "B-spline coefficients must be finite",
        ))
        new(energy_data, mu_data, coefficient_data)
    end
end

"""One linear Legendre-coefficient interpolation interval."""
struct LegendreEnergySpan
    index::Int
    fraction::Float64
    energy_mev::Float64
end

"""
A normalized Legendre angular shape multiplied by a total cross section.

Coefficient rows are stored at increasing centre-of-mass energies in MeV.
Within the tabulated energy domain the coefficients, including ``a0``, are
interpolated linearly. The angular shape is normalized by ``4*pi*a0``. Below
the first energy, the default is an isotropic shape multiplied by the total
cross section; above the final energy the differential cross section is zero.
"""
struct NormalizedLegendreDifferentialCrossSection{T<:TotalCrossSection} <:
       DifferentialCrossSection
    energy_mev::Vector{Float64}
    coefficients::Matrix{Float64}
    total_cross_section::T
    below_minimum_is_isotropic::Bool

    function NormalizedLegendreDifferentialCrossSection(
        energy_mev::AbstractVector{<:Real},
        coefficients::AbstractMatrix{<:Real},
        total_cross_section::T,
        ;
        below_minimum_is_isotropic::Bool=true,
    ) where {T<:TotalCrossSection}
        energies = Float64.(energy_mev)
        coefficient_data = Matrix{Float64}(coefficients)
        length(energies) >= 2 && all(isfinite, energies) &&
            all(diff(energies) .> 0.0) || throw(ArgumentError(
                "Legendre energies must be finite and strictly increasing",
            ))
        size(coefficient_data, 1) == length(energies) &&
            size(coefficient_data, 2) >= 2 &&
            all(isfinite, coefficient_data) || throw(ArgumentError(
                "Legendre coefficients must have one finite row per energy " *
                "and at least a0 and a1",
            ))
        all(coefficient_data[:, 1] .> 0.0) || throw(ArgumentError(
            "Legendre a0 coefficients must be positive",
        ))
        new{T}(
            energies, coefficient_data, total_cross_section,
            below_minimum_is_isotropic,
        )
    end
end

@inline function cubic_bspline_basis(knots::Vector{Float64}, value::Float64)
    coefficient_count = length(knots) - 4
    if value < knots[4] || value > knots[coefficient_count + 1]
        return 0, (0.0, 0.0, 0.0, 0.0)
    end
    span = value == knots[coefficient_count + 1] ? coefficient_count :
        searchsortedlast(knots, value)
    basis = (1.0, 0.0, 0.0, 0.0)
    left = (0.0, 0.0, 0.0, 0.0)
    right = (0.0, 0.0, 0.0, 0.0)
    @inbounds for degree in 1:3
        left = Base.setindex(left, value - knots[span + 1 - degree], degree + 1)
        right = Base.setindex(right, knots[span + degree] - value, degree + 1)
        saved = 0.0
        for basis_index in 0:(degree - 1)
            denominator = right[basis_index + 2] +
                left[degree - basis_index + 1]
            term = denominator == 0.0 ? 0.0 :
                basis[basis_index + 1] / denominator
            basis = Base.setindex(
                basis,
                saved + right[basis_index + 2] * term,
                basis_index + 1,
            )
            saved = left[degree - basis_index + 1] * term
        end
        basis = Base.setindex(basis, saved, degree + 1)
    end
    return span, basis
end

"""Evaluate a tensor B-spline with a precomputed energy basis."""
@inline function differential_cross_section_barns_per_sr(
    model::TensorCubicBSpline, energy_span::Int,
    energy_basis::NTuple{4,Float64}, mu::Float64,
)
    energy_span == 0 && return 0.0
    mu_span, mu_basis = cubic_bspline_basis(model.mu_knots, mu)
    mu_span == 0 && return 0.0
    result = 0.0
    @inbounds for energy_basis_index in 1:4
        coefficient_energy = energy_span + energy_basis_index - 4
        for mu_basis_index in 1:4
            coefficient_mu = mu_span + mu_basis_index - 4
            result += energy_basis[energy_basis_index] *
                mu_basis[mu_basis_index] *
                model.coefficients[coefficient_energy, coefficient_mu]
        end
    end
    # Cubic interpolation can undershoot between non-negative source data.
    # A differential cross section is physical rate data and cannot be negative.
    return isfinite(result) ? max(0.0, result) : result
end

@inline function differential_cross_section_barns_per_sr(
    model::TensorCubicBSpline, energy_mev::Float64, mu::Float64,
)
    energy_span, energy_basis = cubic_bspline_basis(
        model.energy_knots, energy_mev,
    )
    energy_span == 0 && return 0.0
    return differential_cross_section_barns_per_sr(
        model, energy_span, energy_basis, mu,
    )
end

@inline function prepare_differential_cross_section(
    model::TensorCubicBSpline, energy_mev::Float64,
)
    return cubic_bspline_basis(model.energy_knots, energy_mev)
end

@inline function evaluate_differential_cross_section(
    model::TensorCubicBSpline,
    cache::Tuple{Int,NTuple{4,Float64}},
    mu::Float64,
)
    return differential_cross_section_barns_per_sr(
        model, cache[1], cache[2], mu,
    )
end

@inline function prepare_differential_cross_section(
    model::NormalizedLegendreDifferentialCrossSection,
    energy_mev::Float64,
)
    energies = model.energy_mev
    if !isfinite(energy_mev) || energy_mev > energies[end]
        return LegendreEnergySpan(0, 0.0, energy_mev)
    end
    if energy_mev < energies[1]
        return LegendreEnergySpan(-1, 0.0, energy_mev)
    end
    if energy_mev == energies[end]
        return LegendreEnergySpan(length(energies) - 1, 1.0, energy_mev)
    end
    upper = searchsortedfirst(energies, energy_mev)
    if energy_mev == energies[upper]
        return LegendreEnergySpan(upper, 0.0, energy_mev)
    end
    lower = upper - 1
    fraction = (energy_mev - energies[lower]) /
        (energies[upper] - energies[lower])
    return LegendreEnergySpan(lower, fraction, energy_mev)
end

@inline function evaluate_differential_cross_section(
    model::NormalizedLegendreDifferentialCrossSection,
    cache::LegendreEnergySpan,
    mu::Float64,
)
    cache.index == 0 && return 0.0
    -1.0 <= mu <= 1.0 || return 0.0
    sigma_barns = cross_section_barns(
        model.total_cross_section, cache.energy_mev,
    )
    if cache.index == -1
        return model.below_minimum_is_isotropic ?
            max(0.0, sigma_barns / (4.0 * pi)) : 0.0
    end
    lower = cache.index
    upper = lower + 1
    weight = cache.fraction
    coefficients = model.coefficients
    @inbounds begin
        a0 = (1.0 - weight) * coefficients[lower, 1] +
            weight * coefficients[upper, 1]
        a0 > 0.0 || return 0.0
        angular_sum = a0
        p_previous = 1.0
        p_current = mu
        angular_sum += (
            (1.0 - weight) * coefficients[lower, 2] +
            weight * coefficients[upper, 2]
        ) * p_current
        for order in 2:(size(coefficients, 2) - 1)
            p_next = (
                (2order - 1) * mu * p_current -
                (order - 1) * p_previous
            ) / order
            coefficient = (1.0 - weight) * coefficients[lower, order + 1] +
                weight * coefficients[upper, order + 1]
            angular_sum += coefficient * p_next
            p_previous = p_current
            p_current = p_next
        end
    end
    result = sigma_barns * angular_sum / (4.0 * pi * a0)
    return isfinite(result) ? max(0.0, result) : result
end

@inline function differential_cross_section_barns_per_sr(
    model::NormalizedLegendreDifferentialCrossSection,
    energy_mev::Float64,
    mu::Float64,
)
    cache = prepare_differential_cross_section(model, energy_mev)
    return evaluate_differential_cross_section(model, cache, mu)
end
