"""Interface for total cross sections evaluated in barns."""
abstract type TotalCrossSection end

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

"""Tensor-product cubic B-spline used by the DT differential cross section."""
struct TensorCubicBSpline
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
