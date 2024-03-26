using ArgParse
using Printf
using LilGuys



@Base.kwdef mutable struct ParserState
    header = true
    variable = nothing
    data::Dict{String, Vector{Float64}} = Dict()
    sizes::Dict{String, Vector{Int}} = Dict()
    dtypes::Dict{String, String} = Dict()
    finished = false
end


function parse_tsf(filename::String)
    state = ParserState()

    open(filename, "r") do file
        for line in eachline(file)
            parse_line(strip(line), state)
        end
    end

    check_consistency(state)
    return state
end


function parse_line(line, state::ParserState)
    if isempty(line)
        return
    end

    if state.header
        if startswith(line, "set Particles")
            state.header = false
        end
        return
    end

    if startswith(line, "tes")
        state.variable = nothing
        state.finished = true
        return
    end


    elements = split(line)
    if startswith(line, "float")
        state.variable = match(r"(\w+)", elements[2]).captures[1]
        var = state.variable
        state.dtypes[var] = elements[1]
        matches = eachmatch(r"\d+", elements[2])
        sizes = [m.match for m in matches]
        state.sizes[var] = parse.(Int, sizes)
        state.data[var] = Float64[]

        append!(state.data[var], parse.(Float64, elements[3:end]))
        return 
    end

    if state.variable != nothing
        var = state.variable
        append!(state.data[var], parse.(Float64, elements))
    end

end


function check_consistency(state::ParserState)
    if !state.finished
        error("incomplete input")
    end
    for (var, data) in state.data
        expected = prod(state.sizes[var])
        N = length(data)
        if N != expected
            error("incomplete input for $var. expected shape $(expected), got $N")
        end
    end

    for var in ["Mass", "Position", "Velocity"]
        if var ∉ keys(state.data)
            error("no $var variable found")
        end
    end

    masses = state.data["Mass"]
    if ! all(masses .== masses[1])
        error("expected all masses to be the same")
    end
end

function output_to_hdf5(table, filename; verbose=false)
    N = length(table["Mass"])
    mass = table["Mass"][1]
    verbose && println(mass)
end

function inflate(parser::ParserState, var::String)
    return reshape(parser.data[var], reverse(parser.sizes[var])...)
end


function parse_arguments()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "input"
            arg_type=String
            help="input file"
            required=true
        "output"
            arg_type=String
            help="output HDF5 file"
            required=true
    end
    args = parse_args(s)
    return args["input"], args["output"]
end


function main()
    infile, outfile = parse_arguments()
    println("reading ", infile, "\n")

    parsed = parse_tsf(infile)
    N = length(parsed.data["Mass"])
    println("loaded $(N) particles")

    mass = parsed.data["Mass"]
    vel = inflate(parsed, "Velocity")
    pos = inflate(parsed, "Position")

    snap = Snapshot(
        pos,
        vel,
        mass
    )

    println("output: ", outfile, "\n")
    save(outfile, snap) 
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
