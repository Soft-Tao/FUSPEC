using TOML

const BUILTIN_DATA_DIRECTORY = normpath(joinpath(@__DIR__, "..", "data"))

"""A reaction name resolved to particles and both cross-section providers."""
struct ResolvedReaction{T<:TotalCrossSection,D<:DifferentialCrossSection}
    name::String
    reaction::Reaction
    total_cross_section::T
    differential_cross_section::D
end

is_builtin_reaction(name::AbstractString) =
    name == "dd-n3he" || name == "dt-n4he"

function required_toml_field(data::AbstractDict, name::AbstractString,
                             context::AbstractString)
    haskey(data, name) || throw(ArgumentError(
        "$context is missing required field $name",
    ))
    return data[name]
end

@inline function assert_unit_string(data::AbstractDict, name::AbstractString,
                                    expected::AbstractString,
                                    context::AbstractString)
    String(required_toml_field(data, name, context)) == expected ||
        throw(ArgumentError("$context.$name must be $expected"))
    return nothing
end

function toml_float_vector(data::AbstractDict, name::AbstractString,
                           context::AbstractString)
    values = Float64.(required_toml_field(data, name, context))
    all(isfinite, values) || throw(ArgumentError("$context.$name must be finite"))
    return values
end

"""
Load `D(d,n)3He` from its own built-in data file and format.

This entry point may diverge from other built-in reactions when their source
formats or preprocessing requirements change.
"""
function load_dd_n3he_builtin(
    path::AbstractString=joinpath(BUILTIN_DATA_DIRECTORY, "dd-n3he.toml"),
)
    isfile(path) || throw(ArgumentError(
        "built-in reaction data file does not exist: $path",
    ))
    data = TOML.parsefile(path)
    Int(required_toml_field(data, "schema_version", "dd-n3he data")) == 2 ||
        throw(ArgumentError("unsupported dd-n3he data schema_version"))
    String(required_toml_field(data, "reaction_name", "dd-n3he data")) ==
        "dd-n3he" || throw(ArgumentError(
            "dd-n3he data file has the wrong reaction_name",
        ))

    particle_data = required_toml_field(
        data, "particles", "dd-n3he data",
    )
    String.(required_toml_field(
        particle_data, "role_order", "/particles",
    )) == ["r1", "r2", "p3", "p4"] || throw(ArgumentError(
        "dd-n3he /particles.role_order must be [r1, r2, p3, p4]",
    ))
    assert_unit_string(particle_data, "mass_unit", "kg", "/particles")
    atomic_mass_unit = Float64(required_toml_field(
        particle_data, "atomic_mass_unit_kg", "/particles",
    ))
    isfinite(atomic_mass_unit) && atomic_mass_unit > 0.0 ||
        throw(ArgumentError("/particles.atomic_mass_unit_kg must be positive"))
    function dd_particle(role, label, charge_number, mass_number)
        haskey(particle_data, role) || throw(ArgumentError(
            "dd-n3he /particles is missing $role",
        ))
        record = particle_data[role]
        context = "/particles/$role"
        stored_label = String(required_toml_field(record, "label", context))
        mass_u = Float64(required_toml_field(record, "mass_u", context))
        mass_kg = Float64(required_toml_field(record, "mass_kg", context))
        stored_charge = Int(required_toml_field(
            record, "charge_number", context,
        ))
        stored_mass_number = Int(required_toml_field(
            record, "mass_number", context,
        ))
        isfinite(mass_u) && mass_u > 0.0 &&
            isapprox(mass_u * atomic_mass_unit, mass_kg;
                     rtol=1.0e-14, atol=0.0) && stored_charge == charge_number &&
            stored_mass_number == mass_number && stored_label == label ||
            throw(ArgumentError(
                "dd-n3he $context does not match its required particle",
            ))
        return Particle(
            Symbol(stored_label), mass_kg, stored_charge, stored_mass_number,
        )
    end
    r1 = dd_particle("r1", "D", 1, 2)
    r2 = dd_particle("r2", "D", 1, 2)
    p3 = dd_particle("p3", "n", 0, 1)
    p4 = dd_particle("p4", "3He", 2, 3)
    reaction = Reaction(
        :dd_n3he, r1, r2, p3, p4;
        link_to_product1=1, link_to_product2=1,
    )

    total_data = required_toml_field(
        data, "total_cross_section", "dd-n3he data",
    )
    String(required_toml_field(
        total_data, "model", "/total_cross_section",
    )) == "Bosch-Hale" || throw(ArgumentError(
        "dd-n3he total cross section must use Bosch-Hale",
    ))
    assert_unit_string(total_data, "energy_unit", "keV", "/total_cross_section")
    assert_unit_string(
        total_data, "energy_frame", "center-of-mass", "/total_cross_section",
    )
    assert_unit_string(
        total_data, "energy_symbol", "K_CM", "/total_cross_section",
    )
    assert_unit_string(
        total_data, "cross_section_unit", "mb", "/total_cross_section",
    )
    haskey(total_data, "citation") || throw(ArgumentError(
        "dd-n3he total cross section is missing citation",
    ))
    bg = Float64(required_toml_field(
        total_data, "BG_sqrt_keV", "/total_cross_section",
    ))
    branch_data = required_toml_field(
        total_data, "branch", "/total_cross_section",
    )
    length(branch_data) == 1 || throw(ArgumentError(
        "dd-n3he requires exactly one Bosch-Hale branch",
    ))
    branch = branch_data[1]
    total_cross_section = BoschHaleCrossSection([
        BoschHaleBranch(
            Float64(required_toml_field(
                branch, "energy_min_keV", "/total_cross_section/branch[1]",
            )),
            Float64(required_toml_field(
                branch, "energy_max_keV", "/total_cross_section/branch[1]",
            )),
            bg,
            toml_float_vector(
                branch, "A", "/total_cross_section/branch[1]",
            ),
            toml_float_vector(
                branch, "B", "/total_cross_section/branch[1]",
            ),
        ),
    ])

    angular_data = required_toml_field(
        data, "angular_cross_section", "dd-n3he data",
    )
    String(required_toml_field(
        angular_data, "model", "/angular_cross_section",
    )) == "DROSG normalized Legendre expansion" || throw(ArgumentError(
        "dd-n3he angular cross section must use the DROSG expansion",
    ))
    assert_unit_string(
        angular_data, "energy_unit", "MeV", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "energy_frame", "center-of-mass",
        "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "energy_symbol", "K_CM", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "coefficient_unit", "dimensionless",
        "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "energy_interpolation",
        "linear in energy, then normalize by the interpolated a0",
        "/angular_cross_section",
    )
    toml_float_vector(
        angular_data, "mu_range", "/angular_cross_section",
    ) == [-1.0, 1.0] || throw(ArgumentError(
        "dd-n3he /angular_cross_section.mu_range must be [-1, 1]",
    ))
    assert_unit_string(
        angular_data, "below_min_energy", "isotropic",
        "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "above_max_energy", "zero",
        "/angular_cross_section",
    )
    required_toml_field(
        angular_data, "normalized_over_solid_angle", "/angular_cross_section",
    ) === true || throw(ArgumentError(
        "dd-n3he angular coefficients must be normalized",
    ))
    assert_unit_string(
        angular_data, "reference_reactant", "r1", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "reference_product", "p3", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "p4_mu_lookup", "mu -> -mu",
        "/angular_cross_section",
    )
    haskey(angular_data, "citation") || throw(ArgumentError(
        "dd-n3he angular cross section is missing citation",
    ))
    coefficient_count = Int(required_toml_field(
        angular_data, "coefficient_count", "/angular_cross_section",
    ))
    records = required_toml_field(
        angular_data, "records", "/angular_cross_section",
    )
    coefficient_count >= 2 && length(records) >= 2 ||
        throw(ArgumentError("dd-n3he angular data is incomplete"))
    energies = Vector{Float64}(undef, length(records))
    coefficients = Matrix{Float64}(
        undef, length(records), coefficient_count,
    )
    for (index, record) in enumerate(records)
        length(record) == coefficient_count + 1 || throw(DimensionMismatch(
            "dd-n3he angular row $(index - 1) has the wrong length",
        ))
        values = Float64.(record)
        all(isfinite, values) || throw(ArgumentError(
            "dd-n3he angular row $(index - 1) contains non-finite data",
        ))
        energies[index] = values[1]
        coefficients[index, :] .= values[2:end]
    end
    Float64(required_toml_field(
        angular_data, "energy_min_MeV", "/angular_cross_section",
    )) == energies[1] || throw(ArgumentError(
        "dd-n3he angular energy_min_MeV does not match records",
    ))
    Float64(required_toml_field(
        angular_data, "energy_max_MeV", "/angular_cross_section",
    )) == energies[end] || throw(ArgumentError(
        "dd-n3he angular energy_max_MeV does not match records",
    ))
    differential_cross_section = NormalizedLegendreDifferentialCrossSection(
        energies, coefficients, total_cross_section;
        below_minimum_is_isotropic=true,
    )
    return ResolvedReaction(
        "dd-n3he", reaction, total_cross_section,
        differential_cross_section,
    )
end

"""
Load `D(t,n)4He` from its own built-in data file and format.

This entry point may diverge from other built-in reactions when their source
formats or preprocessing requirements change.
"""
function load_dt_n4he_builtin(
    path::AbstractString=joinpath(BUILTIN_DATA_DIRECTORY, "dt-n4he.toml"),
)
    isfile(path) || throw(ArgumentError(
        "built-in reaction data file does not exist: $path",
    ))
    data = TOML.parsefile(path)
    Int(required_toml_field(data, "schema_version", "dt-n4he data")) == 2 ||
        throw(ArgumentError("unsupported dt-n4he data schema_version"))
    String(required_toml_field(data, "reaction_name", "dt-n4he data")) ==
        "dt-n4he" || throw(ArgumentError(
            "dt-n4he data file has the wrong reaction_name",
        ))

    particle_data = required_toml_field(
        data, "particles", "dt-n4he data",
    )
    String.(required_toml_field(
        particle_data, "role_order", "/particles",
    )) == ["r1", "r2", "p3", "p4"] || throw(ArgumentError(
        "dt-n4he /particles.role_order must be [r1, r2, p3, p4]",
    ))
    assert_unit_string(particle_data, "mass_unit", "kg", "/particles")
    atomic_mass_unit = Float64(required_toml_field(
        particle_data, "atomic_mass_unit_kg", "/particles",
    ))
    isfinite(atomic_mass_unit) && atomic_mass_unit > 0.0 ||
        throw(ArgumentError("/particles.atomic_mass_unit_kg must be positive"))
    function dt_particle(role, label, charge_number, mass_number)
        haskey(particle_data, role) || throw(ArgumentError(
            "dt-n4he /particles is missing $role",
        ))
        record = particle_data[role]
        context = "/particles/$role"
        stored_label = String(required_toml_field(record, "label", context))
        mass_u = Float64(required_toml_field(record, "mass_u", context))
        mass_kg = Float64(required_toml_field(record, "mass_kg", context))
        stored_charge = Int(required_toml_field(
            record, "charge_number", context,
        ))
        stored_mass_number = Int(required_toml_field(
            record, "mass_number", context,
        ))
        isfinite(mass_u) && mass_u > 0.0 &&
            isapprox(mass_u * atomic_mass_unit, mass_kg;
                     rtol=1.0e-14, atol=0.0) && stored_charge == charge_number &&
            stored_mass_number == mass_number && stored_label == label ||
            throw(ArgumentError(
                "dt-n4he $context does not match its required particle",
            ))
        return Particle(
            Symbol(stored_label), mass_kg, stored_charge, stored_mass_number,
        )
    end
    r1 = dt_particle("r1", "D", 1, 2)
    r2 = dt_particle("r2", "T", 1, 3)
    p3 = dt_particle("p3", "n", 0, 1)
    p4 = dt_particle("p4", "4He", 2, 4)
    reaction = Reaction(
        :dt_n4he, r1, r2, p3, p4;
        link_to_product1=1, link_to_product2=1,
    )

    total_data = required_toml_field(
        data, "total_cross_section", "dt-n4he data",
    )
    String(required_toml_field(
        total_data, "model", "/total_cross_section",
    )) == "Bosch-Hale" || throw(ArgumentError(
        "dt-n4he total cross section must use Bosch-Hale",
    ))
    assert_unit_string(total_data, "energy_unit", "keV", "/total_cross_section")
    assert_unit_string(
        total_data, "energy_frame", "center-of-mass", "/total_cross_section",
    )
    assert_unit_string(
        total_data, "energy_symbol", "K_CM", "/total_cross_section",
    )
    assert_unit_string(
        total_data, "cross_section_unit", "mb", "/total_cross_section",
    )
    haskey(total_data, "citation") || throw(ArgumentError(
        "dt-n4he total cross section is missing citation",
    ))
    bg = Float64(required_toml_field(
        total_data, "BG_sqrt_keV", "/total_cross_section",
    ))
    branch_data = required_toml_field(
        total_data, "branch", "/total_cross_section",
    )
    length(branch_data) == 2 || throw(ArgumentError(
        "dt-n4he requires two Bosch-Hale branches",
    ))
    branches = BoschHaleBranch[]
    for (index, branch) in enumerate(branch_data)
        context = "/total_cross_section/branch[$index]"
        push!(branches, BoschHaleBranch(
            Float64(required_toml_field(
                branch, "energy_min_keV", context,
            )),
            Float64(required_toml_field(
                branch, "energy_max_keV", context,
            )),
            bg,
            toml_float_vector(branch, "A", context),
            toml_float_vector(branch, "B", context),
        ))
    end
    total_cross_section = BoschHaleCrossSection(branches)

    angular_data = required_toml_field(
        data, "angular_cross_section", "dt-n4he data",
    )
    String(required_toml_field(
        angular_data, "model", "/angular_cross_section",
    )) == "DROSG normalized Legendre expansion" || throw(ArgumentError(
        "dt-n4he angular cross section must use the DROSG expansion",
    ))
    assert_unit_string(
        angular_data, "energy_unit", "MeV", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "energy_frame", "center-of-mass",
        "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "energy_symbol", "K_CM", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "coefficient_unit", "dimensionless",
        "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "energy_interpolation",
        "linear in energy, then normalize by the interpolated a0",
        "/angular_cross_section",
    )
    toml_float_vector(
        angular_data, "mu_range", "/angular_cross_section",
    ) == [-1.0, 1.0] || throw(ArgumentError(
        "dt-n4he /angular_cross_section.mu_range must be [-1, 1]",
    ))
    assert_unit_string(
        angular_data, "below_min_energy", "isotropic",
        "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "above_max_energy", "zero",
        "/angular_cross_section",
    )
    required_toml_field(
        angular_data, "normalized_over_solid_angle", "/angular_cross_section",
    ) === true || throw(ArgumentError(
        "dt-n4he angular coefficients must be normalized",
    ))
    assert_unit_string(
        angular_data, "reference_reactant", "r1", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "reference_product", "p3", "/angular_cross_section",
    )
    assert_unit_string(
        angular_data, "p4_mu_lookup", "mu -> -mu",
        "/angular_cross_section",
    )
    haskey(angular_data, "citation") || throw(ArgumentError(
        "dt-n4he angular cross section is missing citation",
    ))
    coefficient_count = Int(required_toml_field(
        angular_data, "coefficient_count", "/angular_cross_section",
    ))
    records = required_toml_field(
        angular_data, "records", "/angular_cross_section",
    )
    coefficient_count >= 2 && length(records) >= 2 ||
        throw(ArgumentError("dt-n4he angular data is incomplete"))
    energies = Vector{Float64}(undef, length(records))
    coefficients = Matrix{Float64}(
        undef, length(records), coefficient_count,
    )
    for (index, record) in enumerate(records)
        length(record) == coefficient_count + 1 || throw(DimensionMismatch(
            "dt-n4he angular row $(index - 1) has the wrong length",
        ))
        values = Float64.(record)
        all(isfinite, values) || throw(ArgumentError(
            "dt-n4he angular row $(index - 1) contains non-finite data",
        ))
        energies[index] = values[1]
        coefficients[index, :] .= values[2:end]
    end
    Float64(required_toml_field(
        angular_data, "energy_min_MeV", "/angular_cross_section",
    )) == energies[1] || throw(ArgumentError(
        "dt-n4he angular energy_min_MeV does not match records",
    ))
    Float64(required_toml_field(
        angular_data, "energy_max_MeV", "/angular_cross_section",
    )) == energies[end] || throw(ArgumentError(
        "dt-n4he angular energy_max_MeV does not match records",
    ))
    differential_cross_section = NormalizedLegendreDifferentialCrossSection(
        energies, coefficients, total_cross_section;
        below_minimum_is_isotropic=true,
    )
    return ResolvedReaction(
        "dt-n4he", reaction, total_cross_section,
        differential_cross_section,
    )
end

"""
Resolve a supported built-in reaction name through its reaction-specific
loader. The name is matched against a fixed allowlist.
"""
function load_builtin_reaction(
    name::AbstractString;
    data_directory::AbstractString=BUILTIN_DATA_DIRECTORY,
)
    if name == "dd-n3he"
        return load_dd_n3he_builtin(joinpath(data_directory, "dd-n3he.toml"))
    elseif name == "dt-n4he"
        return load_dt_n4he_builtin(joinpath(data_directory, "dt-n4he.toml"))
    end
    throw(ArgumentError("unsupported built-in reaction name $name"))
end
