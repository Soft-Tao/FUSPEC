module FUSPEC

include("constants.jl")
include("particles.jl")
include("kinematics.jl")
include("reactions.jl")
include("rng.jl")
include("distributions.jl")
include("packed_inputs.jl")
include("cross_sections.jl")
include("builtin_reactions.jl")
include("run_config.jl")
include("kinematics_core.jl")
include("reaction_rate.jl")
include("total_spectrum.jl")
include("task_io.jl")
include("capabilities.jl")

export C0, C_SQUARED, ELEMENTARY_CHARGE, MEV
export Particle
export KinematicsMode, Relativistic
export Reaction
export counter_u64, counter_uniform
export Distribution, Maxwellian, BiMaxwellian, TabulatedDistribution,
       PackedDistribution, pack_distribution, sample_momentum
export PackedDistributionBatch, build_reaction_rate_cells
export TotalCrossSection, DifferentialCrossSection, ConstantCrossSection,
       BoschHaleBranch, BoschHaleCrossSection, CubicSplineCrossSection,
       TensorCubicBSpline, NormalizedLegendreDifferentialCrossSection,
       cross_section_barns,
       differential_cross_section_barns_per_sr
export ResolvedReaction, load_builtin_reaction, is_builtin_reaction
export RunConfig
export ReactionRateCell, ReactionRateProblem,
       ReactionRateResult, ReactionRateWorkspace
export calculate_reaction_rate, calculate_reaction_rate!
export SpectrumProblem, EmissionSpectrum, SpectrumWorkspace,
       calculate_total_spectrum, calculate_total_spectrum!,
       normalized_spectrum
export PreparedTask, load_task, write_result
export SharedKinematics, DirectionKinematics, prepare_kinematics,
       direction_kinematics, direction_kinematics_batch!
export kernel_version

end
