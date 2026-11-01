const Vec3 = NTuple{3,Float64}
const ZERO_VEC3 = (0.0, 0.0, 0.0)

@inline dot3(left::Vec3, right::Vec3) =
    left[1] * right[1] + left[2] * right[2] + left[3] * right[3]
@inline norm3(value::Vec3) = sqrt(dot3(value, value))
@inline add3(left::Vec3, right::Vec3) =
    (left[1] + right[1], left[2] + right[2], left[3] + right[3])
@inline subtract3(left::Vec3, right::Vec3) =
    (left[1] - right[1], left[2] - right[2], left[3] - right[3])
@inline scale3(scale::Float64, value::Vec3) =
    (scale * value[1], scale * value[2], scale * value[3])

"""
Quantities shared by every emission direction for one reactant sample.

Both cross sections are evaluated at `center_of_mass_energy`, the kinetic
energy of the relative motion in the centre-of-mass frame, so task files must
tabulate the total and the differential cross section against that same
centre-of-mass kinetic energy.
"""
struct SharedKinematics
    center_of_mass_energy::Float64
    relative_speed::Float64
    feasible::Bool
    selected_mass::Float64
    partner_mass::Float64
    link_to::Int8
    v3_cm::Float64
    v_cm::Float64
    v_cm_vector::Vec3
    v_cm_unit::Vec3
    relativistic_gamma::Float64
    relativistic_g::Float64
    relativistic_m3_energy::Float64
    relativistic_link_momentum::Vec3
    relativistic_link_momentum_norm::Float64
end

"""Direction-dependent product energy, scattering cosine and Jacobian."""
struct DirectionKinematics
    product_energy::Float64
    scattering_cosine::Float64
    jacobian::Float64
    feasible::Bool
end

function product_data(reaction::Reaction, which_product::Integer)
    if which_product == 1
        return reaction.product1.mass, reaction.product2.mass,
               reaction.link_to_product1
    elseif which_product == 2
        return reaction.product2.mass, reaction.product1.mass,
               reaction.link_to_product2
    end
    throw(ArgumentError("which_product must be 1 or 2"))
end

function infeasible_shared(center_of_mass_energy::Float64,
                           relative_speed::Float64,
                           selected_mass::Float64,
                           partner_mass::Float64, link_to::Int8)
    return SharedKinematics(
        center_of_mass_energy, relative_speed, false, selected_mass,
        partner_mass, link_to, 0.0, 0.0, ZERO_VEC3, ZERO_VEC3,
        0.0, 0.0, 0.0, ZERO_VEC3, 0.0,
    )
end

struct FourVector
    energy::Float64
    momentum::Vec3
end

function particle_with_collective(momentum::Vec3, mass::Float64,
                                  collective::Vec3)
    rest_energy = mass * C_SQUARED
    base_energy = sqrt(rest_energy^2 +
                       dot3(momentum, momentum) * C_SQUARED)
    collective_speed_squared = dot3(collective, collective)
    collective_speed_squared < C_SQUARED || throw(ArgumentError(
        "collective velocity must be slower than light",
    ))
    collective_speed_squared == 0.0 && return FourVector(
        base_energy, momentum,
    )
    collective_gamma = 1.0 / sqrt(
        1.0 - collective_speed_squared / C_SQUARED,
    )
    collective_energy = collective_gamma * rest_energy
    collective_momentum = scale3(mass * collective_gamma, collective)
    momentum_dot = dot3(momentum, collective_momentum) * C_SQUARED
    energy = (base_energy * collective_energy + momentum_dot) / rest_energy
    term = momentum_dot / (collective_energy + rest_energy) + base_energy
    boosted_momentum = add3(
        momentum, scale3(term / rest_energy, collective_momentum),
    )
    return FourVector(energy, boosted_momentum)
end

function to_rest_frame(particle::FourVector, frame::FourVector,
                       frame_rest_energy::Float64)
    momentum_dot = dot3(particle.momentum, frame.momentum) * C_SQUARED
    term = momentum_dot / (frame.energy + frame_rest_energy) - particle.energy
    energy = (particle.energy * frame.energy - momentum_dot) /
        frame_rest_energy
    momentum = add3(
        particle.momentum, scale3(term / frame_rest_energy, frame.momentum),
    )
    return FourVector(energy, momentum)
end

function prepare_kinematics(::Relativistic, reaction::Reaction,
                            which_product::Integer, momentum1::Vec3,
                            momentum2::Vec3, collective1::Vec3=ZERO_VEC3,
                            collective2::Vec3=ZERO_VEC3)
    mass1 = reaction.reactant1.mass
    mass2 = reaction.reactant2.mass
    selected_mass, partner_mass, link_to = product_data(
        reaction, which_product,
    )
    particle1 = particle_with_collective(momentum1, mass1, collective1)
    particle2 = particle_with_collective(momentum2, mass2, collective2)
    rest_energy1 = mass1 * C_SQUARED
    rest_energy2 = mass2 * C_SQUARED
    particle21 = to_rest_frame(particle2, particle1, rest_energy1)
    relative_speed = norm3(particle21.momentum) * C_SQUARED /
        particle21.energy

    center_of_mass = FourVector(
        particle1.energy + particle2.energy,
        add3(particle1.momentum, particle2.momentum),
    )
    # Preserve the stage-0 FourMomentum operation order: compute invariant
    # mass first, then convert it back to rest energy. This avoids silently
    # changing the frozen floating-point baseline through cancellation.
    invariant_mass_squared = center_of_mass.energy^2 / C_SQUARED^2 -
        dot3(center_of_mass.momentum, center_of_mass.momentum) / C_SQUARED
    invariant_mass_squared >= 0.0 || throw(ErrorException(
        "negative center-of-mass invariant",
    ))
    invariant_energy = sqrt(invariant_mass_squared) * C_SQUARED
    center_of_mass_energy = invariant_energy - rest_energy1 - rest_energy2
    invariant_energy >= (selected_mass + partner_mass) * C_SQUARED ||
        return infeasible_shared(
            center_of_mass_energy, relative_speed, selected_mass,
            partner_mass, link_to,
        )

    v_cm_vector = scale3(C_SQUARED / center_of_mass.energy,
                         center_of_mass.momentum)
    v_cm = norm3(v_cm_vector)
    v_cm_unit = v_cm > 0.0 ? scale3(1.0 / v_cm, v_cm_vector) : ZERO_VEC3
    selected_rest_energy = selected_mass * C_SQUARED
    partner_rest_energy = partner_mass * C_SQUARED
    momentum_discriminant = invariant_energy^4 + selected_rest_energy^4 +
        partner_rest_energy^4 -
        2.0 * invariant_energy^2 * selected_rest_energy^2 -
        2.0 * invariant_energy^2 * partner_rest_energy^2 -
        2.0 * selected_rest_energy^2 * partner_rest_energy^2
    momentum_discriminant >= 0.0 || throw(ErrorException(
        "negative product momentum discriminant",
    ))
    product_momentum = sqrt(momentum_discriminant) / (2.0 * invariant_energy)
    product_energy3 = (
        invariant_energy^2 + selected_rest_energy^2 - partner_rest_energy^2
    ) / (2.0 * invariant_energy)
    v3_cm = product_momentum / product_energy3 * C0
    relativistic_gamma = center_of_mass.energy / invariant_energy
    relativistic_g = v3_cm > 0.0 ? v_cm / v3_cm : Inf
    linked = link_to == 1 ? particle1 : (link_to == 2 ? particle2 : nothing)
    linked_cm_momentum = ZERO_VEC3
    linked_cm_norm = 0.0
    if linked !== nothing
        linked_cm = to_rest_frame(linked, center_of_mass, invariant_energy)
        linked_cm_momentum = linked_cm.momentum
        linked_cm_norm = norm3(linked_cm_momentum)
    end
    return SharedKinematics(
        center_of_mass_energy, relative_speed, true, selected_mass,
        partner_mass, link_to, v3_cm, v_cm, v_cm_vector, v_cm_unit,
        relativistic_gamma, relativistic_g,
        selected_rest_energy, linked_cm_momentum, linked_cm_norm,
    )
end

function direction_kinematics(::Relativistic, shared::SharedKinematics,
                              emission_direction::Vec3)
    shared.feasible || return DirectionKinematics(0.0, 0.0, 0.0, false)
    cos_phi = clamp(dot3(emission_direction, shared.v_cm_unit), -1.0, 1.0)
    sin_phi_squared = max(0.0, 1.0 - cos_phi^2)
    gamma_squared = shared.relativistic_gamma^2
    temp2 = gamma_squared * (1.0 - shared.relativistic_g^2) *
        sin_phi_squared + cos_phi^2
    temp2 >= 0.0 || return DirectionKinematics(0.0, 0.0, 0.0, false)
    cos_theta1 = (
        cos_phi * sqrt(max(0.0, temp2)) - gamma_squared *
        shared.relativistic_g * sin_phi_squared
    ) / (gamma_squared * sin_phi_squared + cos_phi^2)
    v3_cm_x = shared.v3_cm * cos_theta1
    v3_cm_y = shared.v3_cm * sqrt(max(0.0, 1.0 - cos_theta1^2))
    numerator = sqrt(
        gamma_squared * (v3_cm_x + shared.v_cm)^2 + v3_cm_y^2,
    )
    denominator = shared.relativistic_gamma *
        (1.0 + v3_cm_x * shared.v_cm / C_SQUARED)
    product_lab_speed = numerator / denominator
    product_gamma = 1.0 / sqrt(max(
        1.0 - product_lab_speed^2 / C_SQUARED, 1e-15,
    ))
    product_energy = (product_gamma - 1.0) *
        shared.relativistic_m3_energy
    a = gamma_squared * (shared.relativistic_g + cos_theta1)^2 +
        1.0 - cos_theta1^2
    jacobian = a^1.5 / (
        shared.relativistic_gamma *
        (1.0 + shared.relativistic_g * cos_theta1)
    )

    total_energy = product_energy + shared.relativistic_m3_energy
    product_momentum_magnitude = sqrt(max(
        total_energy^2 / C_SQUARED -
        shared.selected_mass^2 * C_SQUARED, 0.0,
    ))
    product_momentum = scale3(product_momentum_magnitude, emission_direction)
    beta = scale3(1.0 / C0, shared.v_cm_vector)
    beta_squared = dot3(beta, beta)
    product_cm_momentum = product_momentum
    if beta_squared > 1e-12
        factor = (shared.relativistic_gamma - 1.0) / beta_squared
        beta_dot_momentum = dot3(product_momentum, beta)
        product_cm_momentum = subtract3(
            add3(product_momentum,
                 scale3(factor * beta_dot_momentum, beta)),
            scale3(shared.relativistic_gamma * total_energy / C0, beta),
        )
    end
    product_cm_norm = norm3(product_cm_momentum)
    scattering_cosine = shared.relativistic_link_momentum_norm > 0.0 ?
        clamp(dot3(product_cm_momentum, shared.relativistic_link_momentum) /
              (product_cm_norm * shared.relativistic_link_momentum_norm),
              -1.0, 1.0) : 0.0
    return DirectionKinematics(
        product_energy, scattering_cosine, jacobian, true,
    )
end

"""
Batch direction kinematics for all emission directions of one sample.

Computes the same quantities as [`direction_kinematics`](@ref) for every
direction in `emission_directions` and writes them into caller-provided,
preallocated arrays.  Each direction is independent, so the loop has no
cross-lane reduction; the per-direction floating-point operation order is
identical to the scalar path.  A direction whose `temp2` discriminant is
negative is marked `feasible[i] == false` and its energy/cosine/Jacobian
entries are unspecified (the caller must skip them via `feasible`).
"""
function direction_kinematics_batch!(
    product_energy::AbstractVector{Float64},
    scattering_cosine::AbstractVector{Float64},
    jacobian::AbstractVector{Float64},
    feasible::AbstractVector{Bool},
    ::Relativistic,
    shared::SharedKinematics,
    emission_directions::Vector{Vec3},
)
    directions = length(emission_directions)
    @inbounds if !shared.feasible
        for direction_index in 1:directions
            product_energy[direction_index] = 0.0
            scattering_cosine[direction_index] = 0.0
            jacobian[direction_index] = 0.0
            feasible[direction_index] = false
        end
        return nothing
    end

    gamma_squared = shared.relativistic_gamma^2
    relativistic_g = shared.relativistic_g
    v3_cm = shared.v3_cm
    v_cm = shared.v_cm
    gamma = shared.relativistic_gamma
    m3_energy = shared.relativistic_m3_energy
    selected_mass = shared.selected_mass
    v_cm_unit = shared.v_cm_unit
    beta = scale3(1.0 / C0, shared.v_cm_vector)
    beta_squared = dot3(beta, beta)
    link_momentum = shared.relativistic_link_momentum
    link_norm = shared.relativistic_link_momentum_norm

    @inbounds @simd for direction_index in 1:directions
        dir = emission_directions[direction_index]
        cos_phi = clamp(dot3(dir, v_cm_unit), -1.0, 1.0)
        sin_phi_squared = max(0.0, 1.0 - cos_phi^2)
        temp2 = gamma_squared * (1.0 - relativistic_g^2) *
            sin_phi_squared + cos_phi^2
        feasible[direction_index] = temp2 >= 0.0
        cos_theta1 = (
            cos_phi * sqrt(max(0.0, temp2)) - gamma_squared *
            relativistic_g * sin_phi_squared
        ) / (gamma_squared * sin_phi_squared + cos_phi^2)
        v3_cm_x = v3_cm * cos_theta1
        v3_cm_y = v3_cm * sqrt(max(0.0, 1.0 - cos_theta1^2))
        numerator = sqrt(
            gamma_squared * (v3_cm_x + v_cm)^2 + v3_cm_y^2,
        )
        denominator = gamma * (1.0 + v3_cm_x * v_cm / C_SQUARED)
        product_lab_speed = numerator / denominator
        product_gamma = 1.0 / sqrt(max(
            1.0 - product_lab_speed^2 / C_SQUARED, 1e-15,
        ))
        product_energy_value = (product_gamma - 1.0) * m3_energy
        a = gamma_squared * (relativistic_g + cos_theta1)^2 +
            1.0 - cos_theta1^2
        jacobian[direction_index] = a^1.5 / (
            gamma * (1.0 + relativistic_g * cos_theta1)
        )
        product_energy[direction_index] = product_energy_value

        total_energy = product_energy_value + m3_energy
        product_momentum_magnitude = sqrt(max(
            total_energy^2 / C_SQUARED -
            selected_mass^2 * C_SQUARED, 0.0,
        ))
        product_momentum = scale3(product_momentum_magnitude, dir)
        product_cm_momentum = product_momentum
        if beta_squared > 1e-12
            factor = (gamma - 1.0) / beta_squared
            beta_dot_momentum = dot3(product_momentum, beta)
            product_cm_momentum = subtract3(
                add3(product_momentum,
                     scale3(factor * beta_dot_momentum, beta)),
                scale3(gamma * total_energy / C0, beta),
            )
        end
        product_cm_norm = norm3(product_cm_momentum)
        scattering_cosine[direction_index] = link_norm > 0.0 ?
            clamp(dot3(product_cm_momentum, link_momentum) /
                  (product_cm_norm * link_norm), -1.0, 1.0) : 0.0
    end
    return nothing
end
