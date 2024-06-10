#!/usr/bin/env julia

using ArgParse
using Printf
using LilGuys
using DataFrames, CSV

include("tsf_parser.jl")


function check_is_profile(state::ParserState)
    for var in ["Radius", "Mass", "Density"]
        if var ∉ keys(state.data)
            error("no $var variable found")
        end
    end

end


function output_to_csv(table, filename; verbose=false)
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
            help="output csv file"
            required=false
            default=nothing
    end
    args = parse_args(s)

    if args["output"] == nothing
        args["output"] = replace(args["input"], r"\.[a-z]+$" => ".csv")
    end

    return args["input"], args["output"]
end


function main()
    infile, outfile = parse_arguments()
    println("reading ", infile, "\n")

    parsed = parse_tsf(infile, set="GeneralSphericalProfile")
    check_is_profile(parsed)

    df = DataFrame(
        Radius = parsed.data["Radius"],
        Mass = parsed.data["Mass"],
        Density = parsed.data["Density"]
    )

    println(length(df.Radius), 
        " points")

    CSV.write(outfile, df)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
