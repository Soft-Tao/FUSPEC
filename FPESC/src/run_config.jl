"""Validated configuration shared by FPESC Monte Carlo kernels."""
struct RunConfig{K<:KinematicsMode}
    samples::UInt64
    sample_offset::UInt64
    seed::UInt64
    n_workers::Int
    kinematics::K
    collective_velocity1::NTuple{3,Float64}
    collective_velocity2::NTuple{3,Float64}
end

function velocity_tuple(name::AbstractString, value)
    length(value) == 3 || throw(ArgumentError("$name must have length 3"))
    velocity = (Float64(value[1]), Float64(value[2]), Float64(value[3]))
    all(isfinite, velocity) || throw(ArgumentError("$name must be finite"))
    sum(component^2 for component in velocity) < C0^2 ||
        throw(ArgumentError("$name must be slower than light"))
    return velocity
end

function RunConfig(; samples::Integer=100_000,
                   sample_offset::Integer=0,
                   seed::Integer=0,
                   n_workers::Integer=1,
                   kinematics::KinematicsMode=Relativistic(),
                   collective_velocity1=(0.0, 0.0, 0.0),
                   collective_velocity2=(0.0, 0.0, 0.0))
    samples >= 0 || throw(ArgumentError("samples must be non-negative"))
    sample_offset >= 0 || throw(ArgumentError(
        "sample_offset must be non-negative",
    ))
    seed >= 0 || throw(ArgumentError("seed must be non-negative"))
    n_workers >= 1 || throw(ArgumentError("n_workers must be at least 1"))
    samples <= typemax(UInt64) || throw(ArgumentError("samples exceed UInt64"))
    sample_offset <= typemax(UInt64) || throw(ArgumentError(
        "sample_offset exceeds UInt64",
    ))
    seed <= typemax(UInt64) || throw(ArgumentError("seed exceeds UInt64"))
    samples <= typemax(UInt64) - sample_offset || throw(ArgumentError(
        "samples + sample_offset overflow UInt64",
    ))
    n_workers <= Threads.nthreads() || throw(ArgumentError(
        "n_workers=$(n_workers) exceeds Julia thread capacity " *
        "$(Threads.nthreads())",
    ))
    velocity1 = velocity_tuple("collective_velocity1", collective_velocity1)
    velocity2 = velocity_tuple("collective_velocity2", collective_velocity2)
    return RunConfig(UInt64(samples), UInt64(sample_offset), UInt64(seed),
                     Int(n_workers), kinematics, velocity1, velocity2)
end
