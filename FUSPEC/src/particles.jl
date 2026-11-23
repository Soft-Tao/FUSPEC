"""A particle species expressed in the SI units used by FUSPEC."""
struct Particle
    name::Symbol
    mass::Float64
    charge_number::Int8
    mass_number::Int16

    function Particle(name::Symbol, mass::Real, charge_number::Integer,
                      mass_number::Integer=0)
        isfinite(mass) && mass > 0 || throw(ArgumentError(
            "particle mass must be finite and positive",
        ))
        typemin(Int8) <= charge_number <= typemax(Int8) || throw(ArgumentError(
            "charge number does not fit Int8",
        ))
        0 <= mass_number <= typemax(Int16) || throw(ArgumentError(
            "mass number must be non-negative and fit Int16",
        ))
        new(name, Float64(mass), Int8(charge_number), Int16(mass_number))
    end
end

Base.:(==)(left::Particle, right::Particle) =
    left.name == right.name && left.mass == right.mass &&
    left.charge_number == right.charge_number &&
    left.mass_number == right.mass_number
