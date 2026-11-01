"""
Regression tests for the FPESC physics kernels.

The cases are small enough to run in a few seconds, but they pin the two
properties that are easy to break silently and expensive to notice:

* the emission directions must be expressed in the same field-aligned frame as
  the sampled reactant momenta, because `sample_momentum` uses the local
  magnetic field as the polar axis (TaskPackageProtocol section 7), and
* the Monte Carlo result must not depend on the worker count.

Run them with `julia -t 4 --project=FPESC -e 'using Pkg; Pkg.test()'`; the
worker-count test is skipped when Julia was started with a single thread.
"""

using Test
using FPESC

const AMU = 1.66053906660e-27
const MU_B = [-1.0, -0.5, 0.0, 0.5, 1.0]

mev(value) = value * 1.0e6 * ELEMENTARY_CHARGE

function dt_reaction()
    deuteron = Particle(:d, 2.013553198 * AMU, 1, 2)
    triton = Particle(:t, 3.015500688 * AMU, 1, 3)
    neutron = Particle(:n, 1.008664923 * AMU, 0, 1)
    alpha = Particle(:a, 4.001506090 * AMU, 2, 4)
    return Reaction(
        :dt, deuteron, triton, neutron, alpha;
        link_to_product1=1, link_to_product2=1,
    )
end

"""
Forward-peaked differential cross section on a 12 MeV x [-1, 1] domain.

The coefficient ramp is linear in the mu index, so the spline reproduces
`sigma(mu) ~ 1 + 0.5*mu` and the forward direction carries the largest weight.
"""
function anisotropic_tensor()
    energy_knots = [0.0, 0.0, 0.0, 0.0, 2.0, 4.0, 6.0, 8.0, 10.0,
                    12.0, 12.0, 12.0, 12.0]
    mu_knots = [-1.0, -1.0, -1.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
    rows = length(energy_knots) - 4
    columns = length(mu_knots) - 4
    coefficients = [
        1.0 + 0.5 * (-1.0 + 2.0 * (column - 1) / (columns - 1))
        for _ in 1:rows, column in 1:columns
    ]
    return TensorCubicBSpline(energy_knots, mu_knots, coefficients)
end

"""One cell of 1 MeV deuterium streaming along the frame axis into cold tritium."""
function beam_like_cell()
    energies = mev.([0.5, 1.0, 1.5, 2.0])
    values = repeat([0.02 0.2 1.0], 4, 1)
    return ReactionRateCell(
        1.0e19, 1.0e19, 1.0,
        pack_distribution(
            TabulatedDistribution(energies, [-1.0, 0.0, 1.0], values),
        ),
        pack_distribution(Maxwellian(mev(0.02))),
    )
end

function beam_like_problem()
    cell = beam_like_cell()
    sampling = ReactionRateProblem(
        dt_reaction(), [cell], ConstantCrossSection(0.0);
        energy_range1=(mev(0.5), mev(2.0)), energy_range2=(0.0, mev(0.2)),
        mu_range1=(-1.0, 1.0), mu_range2=(-1.0, 1.0), which_product=1,
    )
    return SpectrumProblem(
        sampling, anisotropic_tensor();
        emission_mus=MU_B, energy_range=(mev(10.0), mev(20.0)), bins=8,
    )
end

"""
`f(E) = 2 + 3E - E^2 + 0.5E^3` written as a two-interval SciPy power basis.

Within an interval the coefficients are (cubic, quadratic, linear, constant)
at `dx = E - x_i`, so they are `(0.5, -1 + 1.5x, 3 - 2x + 1.5x^2, f(x))`.
"""
function analytic_total_cross_section()
    return CubicSplineCrossSection(
        [0.0, 1.0, 2.0], [0.5 0.5; -1.0 0.5; 3.0 2.5; 2.0 4.5],
    )
end

analytic_cubic(energy) = 2.0 + 3.0 * energy - energy^2 + 0.5 * energy^3

"""Reaction-rate counterpart of [`beam_like_problem`](@ref)."""
function beam_like_rate_problem()
    return ReactionRateProblem(
        dt_reaction(), [beam_like_cell()], analytic_total_cross_section();
        energy_range1=(mev(0.5), mev(2.0)), energy_range2=(0.0, mev(0.2)),
        mu_range1=(-1.0, 1.0), mu_range2=(-1.0, 1.0), which_product=1,
    )
end

"""
Isotropic plasma, built twice.

`Maxwellian` pairs take the isotropic shortcut in the kernel; an equivalent
`BiMaxwellian` pair with parallel == perpendicular temperature does not, so the
two problems exercise the two code paths on the same physics.
"""
function isotropic_problem(parallel_temperature, perpendicular_temperature)
    temperature = mev(0.01)
    first = parallel_temperature === nothing ?
        Maxwellian(temperature) :
        BiMaxwellian(parallel_temperature, perpendicular_temperature)
    second = Maxwellian(temperature)
    cell = ReactionRateCell(
        1.0e19, 1.0e19, 1.0, pack_distribution(first), pack_distribution(second),
    )
    sampling = ReactionRateProblem(
        dt_reaction(), [cell], ConstantCrossSection(0.0);
        energy_range1=(0.0, mev(0.2)), energy_range2=(0.0, mev(0.2)),
        mu_range1=(-1.0, 1.0), mu_range2=(-1.0, 1.0), which_product=1,
    )
    return SpectrumProblem(
        sampling, anisotropic_tensor();
        emission_mus=MU_B, energy_range=(mev(13.0), mev(15.0)), bins=8,
    )
end

"""
Frozen accumulator of [`beam_like_problem`](@ref) for 200_000 samples, seed
12345. Any change to the sampler, the kinematics, the cross-section lookup or
the accumulation order shows up here. Regenerate it only together with a
deliberate physics change.
"""
const BEAM_GOLDEN = [
    1.5068907620972094e43 2.518446428016522e46 1.1880793861657409e47 6.979091198784093e46 3.3697882948604606e46 1.322342131027946e46 1.453768881272336e45 6.68457713493535e39;
    6.53456655014754e42 1.5392814333460148e46 8.8891349559189e46 9.488856841072701e46 7.015410860425703e46 3.225664069733341e46 5.194993061191817e45 8.570671179280521e41;
    3.783640277221294e42 8.113939245390406e45 6.245432230976151e46 9.92876472644522e46 1.063235646436249e47 6.54372605484466e46 1.1649037598279196e46 6.251560969017792e41;
    1.1231371055320928e42 3.598329441772941e45 3.726423603519011e46 8.91746170596064e46 1.3745367330860746e47 1.0944567123326619e47 2.508135231349343e46 2.7075121794864007e42;
    2.732872940782582e41 1.2168493966153427e45 1.7127265808087597e46 5.881274602541087e46 1.6599997768284488e47 1.6884209650021914e47 4.1102394445596066e46 4.926527836727582e42
]

"""Beam-thermal samples, identical for every worker count."""
function run_beam_problem(n_workers)
    problem = beam_like_problem()
    config = RunConfig(
        samples=200_000, seed=12345, n_workers=n_workers,
        kinematics=Relativistic(),
    )
    return calculate_total_spectrum!(EmissionSpectrum(problem), problem, config)
end

"""Reaction-rate counterpart of [`run_beam_problem`](@ref)."""
function run_beam_rate_problem(n_workers)
    problem = beam_like_rate_problem()
    config = RunConfig(
        samples=500_000, seed=4242, n_workers=n_workers,
        kinematics=Relativistic(),
    )
    return calculate_reaction_rate!(ReactionRateResult(1), problem, config)
end

"""Frozen reaction rate and sample weight for [`beam_like_rate_problem`](@ref)."""
const RATE_GOLDEN = 4.402366809013499e-21
const RATE_WEIGHT_GOLDEN = 8.18184830154786e30

@testset "FPESC" begin
    @testset "counter RNG is deterministic and unbiased in its dimensions" begin
        @test counter_u64(UInt64(7), UInt64(3), UInt32(1)) ==
            counter_u64(UInt64(7), UInt64(3), UInt32(1))
        @test counter_u64(UInt64(7), UInt64(3), UInt32(1)) !=
            counter_u64(UInt64(7), UInt64(3), UInt32(2))
        @test counter_u64(UInt64(7), UInt64(3), UInt32(1)) !=
            counter_u64(UInt64(7), UInt64(4), UInt32(1))
        draws = [counter_uniform(UInt64(20260922), UInt64(index), UInt32(1))
                 for index in 0:9999]
        @test all(0.0 .<= draws .< 1.0)
        @test 0.49 < sum(draws) / length(draws) < 0.51
    end

    @testset "relative speed and centre-of-mass energy" begin
        reaction = dt_reaction()
        first_energy = mev(0.1)
        second_energy = mev(0.1)
        mass1 = reaction.reactant1.mass
        mass2 = reaction.reactant2.mass
        momentum1 = sqrt(2.0 * mass1 * first_energy)
        momentum2 = sqrt(2.0 * mass2 * second_energy)
        shared = prepare_kinematics(
            Relativistic(), reaction, 1,
            (-momentum1, 0.0, 0.0), (momentum2, 0.0, 0.0),
        )

        # Head-on collinear motion: the relative speed follows the relativistic
        # velocity-addition rule, which pins the unit conventions of the
        # stored momenta (m*v, not p*c).
        speed1 = momentum1 * C_SQUARED /
            sqrt((mass1 * C_SQUARED)^2 + (momentum1 * C0)^2)
        speed2 = momentum2 * C_SQUARED /
            sqrt((mass2 * C_SQUARED)^2 + (momentum2 * C0)^2)
        @test shared.relative_speed ≈
            (speed1 + speed2) / (1.0 + speed1 * speed2 / C_SQUARED) rtol=1e-12

        rest_energy1 = mass1 * C_SQUARED
        rest_energy2 = mass2 * C_SQUARED
        energy1 = sqrt(rest_energy1^2 + (momentum1 * C0)^2)
        energy2 = sqrt(rest_energy2^2 + (momentum2 * C0)^2)
        invariant = (energy1 + energy2)^2 -
            (momentum2 - momentum1)^2 * C_SQUARED
        # The kernel forms this energy by subtracting two ~1e-9 J rest energies,
        # so only a handful of digits of the ~1e-14 J result survive.
        @test shared.center_of_mass_energy ≈
            sqrt(invariant) - rest_energy1 - rest_energy2 rtol=1e-9
    end

    @testset "emission directions share the sampler's polar axis" begin
        problem = beam_like_problem()
        @test length(problem.emission_directions) == length(MU_B)
        for (index, mu) in enumerate(MU_B)
            direction = problem.emission_directions[index]
            @test sum(abs2, direction) ≈ 1.0 rtol=1e-15
            # The emission direction must resolve the same polar axis that
            # `sample_momentum` uses, otherwise the kernel would mix the
            # field-aligned frame with the laboratory basis.
            @test direction[3] ≈ mu rtol=1e-15

            momentum, _ = sample_momentum(
                pack_distribution(Maxwellian(mev(0.02))),
                dt_reaction().reactant1.mass,
                (0.0, mev(0.2)), (mu, mu), UInt64(12345), UInt64(0), UInt32(1),
            )
            @test momentum[3] / sqrt(sum(abs2, momentum)) ≈ mu rtol=1e-14
        end
    end

    @testset "tensor B-spline evaluation" begin
        energy_knots = [0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
                        6.0, 6.0, 6.0, 6.0]
        mu_knots = [-1.0, -1.0, -1.0, -1.0, -0.5, 0.0, 0.5, 1.0, 1.0, 1.0, 1.0]
        rows = length(energy_knots) - 4
        columns = length(mu_knots) - 4
        constant = TensorCubicBSpline(energy_knots, mu_knots, fill(2.5, rows, columns))
        # Partition of unity: constant coefficients reproduce the constant.
        for energy in (0.0, 0.4, 1.5, 3.7, 6.0), mu in (-1.0, -0.3, 0.0, 0.8, 1.0)
            @test differential_cross_section_barns_per_sr(constant, energy, mu) ≈
                2.5 rtol=1e-12
        end
        # Outside the knot domain the tensor returns exactly zero.
        @test differential_cross_section_barns_per_sr(constant, -0.1, 0.5) == 0.0
        @test differential_cross_section_barns_per_sr(constant, 6.1, 0.5) == 0.0
        @test differential_cross_section_barns_per_sr(constant, 3.0, 1.5) == 0.0
        # Cubic undershoot between non-negative samples is clamped away.
        # Coefficients of 1.0 reproduce a flat 1.0, and the single large
        # negative coefficient carves a dip below zero near (E, mu) = (2, -1).
        dipping = fill(1.0, rows, columns)
        dipping[3, 1] = -6.0
        undershoot = TensorCubicBSpline(
            energy_knots, mu_knots, dipping,
        )
        sampled = [
            differential_cross_section_barns_per_sr(undershoot, energy, mu)
            for energy in 0.0:0.05:6.0, mu in -1.0:0.05:1.0
        ]
        @test all(sampled .>= 0.0)
        # The negative lobe is clamped to exactly zero, not merely absent.
        @test minimum(sampled) == 0.0
        @test maximum(sampled) > 0.0
    end

    @testset "total cross-section spline and its extrapolation contract" begin
        cross_section = analytic_total_cross_section()
        for energy in (0.0, 0.25, 0.5, 1.0, 1.4, 2.0)
            @test cross_section_barns(cross_section, energy) ≈
                analytic_cubic(energy) rtol=1e-14
        end
        # Schema v1 continues the first/last polynomial outside the knots
        # instead of returning zero.
        for energy in (-0.5, -0.01, 2.01, 3.0)
            @test cross_section_barns(cross_section, energy) ≈
                analytic_cubic(energy) rtol=1e-12
        end
    end

    @testset "reaction-rate golden" begin
        result = run_beam_rate_problem(1)
        @test result.samples_completed == UInt64(500_000)
        @test result.infeasible_samples == UInt64(0)
        @test result.weight_sums[1] ≈ RATE_WEIGHT_GOLDEN rtol=1e-12
        @test result.reaction_rate[1] ≈ RATE_GOLDEN rtol=1e-12
    end

    @testset "isotropic shortcut matches the explicit direction loop" begin
        shortcut = isotropic_problem(nothing, nothing)
        explicit = isotropic_problem(mev(0.01), mev(0.01))
        @test shortcut.is_isotropic
        @test !explicit.is_isotropic

        config = RunConfig(
            samples=1_000_000, seed=20260922, n_workers=1, kinematics=Relativistic(),
        )
        representative = shortcut.isotropic_direction_index
        shortcut_result = calculate_total_spectrum!(
            EmissionSpectrum(shortcut), shortcut, config,
        )
        # The shortcut computes one direction and copies it, exactly.
        for index in axes(shortcut_result.accumulator, 2)
            @test shortcut_result.accumulator[1, index, :] ==
                shortcut_result.accumulator[1, representative, :]
        end
        @test any(shortcut_result.accumulator[1, representative, :] .> 0.0)

        # An isotropic plasma produces a direction-independent spectrum, so the
        # shortcut must agree with the full per-direction loop in expectation.
        # The comparison is statistical: the thermal estimator has a small
        # effective sample size, so 1e6 samples still scatter by ~5% between
        # directions (measured over several seeds). The tolerance below is a
        # smoke test that would still catch a frame or direction error, which
        # shifts these spectra by tens of percent.
        explicit_result = calculate_total_spectrum!(
            EmissionSpectrum(explicit), explicit,
            RunConfig(samples=1_000_000, seed=20260922, n_workers=1,
                      kinematics=Relativistic()),
        )
        for index in axes(explicit_result.accumulator, 2)
            @test isapprox(
                explicit_result.accumulator[1, index, :] ./
                    explicit_result.weight_sums[1],
                shortcut_result.accumulator[1, representative, :] ./
                    shortcut_result.weight_sums[1];
                rtol=0.25,
            )
        end
    end

    @testset "beam-like spectrum: Doppler ordering and frozen values" begin
        problem = beam_like_problem()
        result = run_beam_problem(1)
        accumulator = result.accumulator[1, :, :]
        centres = [
            problem.energy_min + (bin - 0.5) * problem.bin_width
            for bin in axes(accumulator, 2)
        ]

        means = [sum(centres .* accumulator[index, :]) /
                 sum(accumulator[index, :]) for index in eachindex(MU_B)]
        totals = vec(sum(accumulator, dims=2))
        # Forward emission must be harder and stronger than backward emission.
        @test issorted(means)
        @test issorted(totals)
        @test means[end] > 1.05 * means[1]

        @test all(isapprox.(accumulator, BEAM_GOLDEN; rtol=1e-10))
    end

    @testset "worker-count invariance" begin
        if Threads.nthreads() == 1
            @info "worker-count invariance is skipped: run with julia -t N, N >= 2"
        else
            reference = run_beam_problem(1).accumulator
            rate_reference = run_beam_rate_problem(1).reaction_rate
            for workers in (min(2, Threads.nthreads()), Threads.nthreads())
                @test run_beam_problem(workers).accumulator == reference
                @test run_beam_rate_problem(workers).reaction_rate == rate_reference
            end
        end
    end
end
