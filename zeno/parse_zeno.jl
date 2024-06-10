using ArgParse
using Printf
using LilGuys


include("tsf_parser.jl")


function check_is_snapshot(state::ParserState)
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
    check_is_snapshot(parsed)

    N = length(parsed.data["Mass"])
    println("loaded $(N) particles")

    mass = parsed.data["Mass"]
    if all(mass .== mass[1])
        mass = mass[1]
    end
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
