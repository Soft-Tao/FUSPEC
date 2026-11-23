#!/usr/bin/env julia

using FUSPEC

struct CLIOptions
    input::String
    output::String
    samples::UInt64
end

function print_usage(io::IO=stdout)
    println(io, "Usage:")
    println(io, "  julia -t N --project=FUSPEC examples/scripts/RunTask.jl \\")
    println(io, "    --input TASK.h5 --output OUTPUT.h5 --samples COUNT")
    println(io)
    println(io, "COUNT is the total number of Monte Carlo samples for this run.")
    println(io, "The worker count is Threads.nthreads(), selected by julia -t N.")
end

function parse_samples(value::AbstractString)
    samples = try
        parse(UInt64, value)
    catch
        throw(ArgumentError("--samples must be an integer in [1, 2^64-1]"))
    end
    samples > 0 || throw(ArgumentError(
        "--samples must be an integer in [1, 2^64-1]",
    ))
    return samples
end

function parse_cli(arguments::Vector{String})
    any(argument -> argument in ("-h", "--help"), arguments) && return nothing
    values = Dict{String,String}()
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        option in ("--input", "--output", "--samples") ||
            throw(ArgumentError("unknown CLI option: $option"))
        haskey(values, option) && throw(ArgumentError(
            "CLI option $option was provided more than once",
        ))
        index < length(arguments) || throw(ArgumentError(
            "CLI option $option requires a value",
        ))
        values[option] = arguments[index + 1]
        index += 2
    end
    for option in ("--input", "--output", "--samples")
        haskey(values, option) || throw(ArgumentError(
            "missing required CLI option $option",
        ))
    end
    input = abspath(values["--input"])
    output = abspath(values["--output"])
    input == output && throw(ArgumentError(
        "--input and --output must refer to different files",
    ))
    return CLIOptions(input, output, parse_samples(values["--samples"]))
end

function atomic_write_result(output::AbstractString, task, result, config)
    parent = dirname(output)
    isdir(parent) || mkpath(parent)
    temporary, stream = mktemp(parent; cleanup=false)
    close(stream)
    try
        write_result(temporary, task, result, config)
        mv(temporary, output; force=true)
    finally
        isfile(temporary) && rm(temporary)
    end
    return output
end

function run(options::CLIOptions)
    println("Loading FUSPEC task: ", options.input)
    task = load_task(options.input)
    workers = Threads.nthreads()
    config = RunConfig(
        samples=options.samples,
        seed=task.seed,
        n_workers=workers,
        kinematics=Relativistic(),
    )
    println("Calculation kind: ", replace(String(task.calculation_kind), "_" => "-"))
    println("Samples: ", options.samples)
    println("Julia threads/workers: ", workers)
    result = nothing
    elapsed = @elapsed begin
        result = if task.calculation_kind == :product_emission_spectrum
            calculate_total_spectrum(task.problem, config)
        elseif task.calculation_kind == :reaction_rate_profile
            calculate_reaction_rate(task.problem, config)
        else
            error("unsupported loaded calculation kind $(task.calculation_kind)")
        end
    end
    println("Kernel elapsed seconds: ", round(elapsed; digits=3))
    atomic_write_result(options.output, task, result, config)
    println("Result written: ", options.output)
    return nothing
end

function main(arguments::Vector{String}=ARGS)
    options = try
        parse_cli(arguments)
    catch error
        println(stderr, "Argument error: ", sprint(showerror, error))
        print_usage(stderr)
        return 2
    end
    if options === nothing
        print_usage()
        return 0
    end
    try
        run(options)
        return 0
    catch error
        println(stderr, "FUSPEC task failed: ", sprint(showerror, error))
        Base.display_error(Base.catch_stack())
        return 1
    end
end

exit(main())
