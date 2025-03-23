#!/bin/env julia

using ArgParse
using LilGuys
import PythonCall # for fits
import TOML
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")


function get_args()

    s = ArgParseSettings(description="applies the given filters to a data file and writes the selection")


    @add_arg_table s begin
        "input"
            help = "path to TOML file for filter & profile settings"
            required = true
        "--output", "-o"
            help = "output file"
    end

    args = parse_args(s)

    if isnothing(args["output"])
        filebase, _ = splitext(args["input"])
        args["output"] = "$(filebase)_sample.fits"
    end

    return args
end

function main()
    args = get_args()

    params = read_paramfile(args["input"])

    profile_kwargs = pop!(params, "profile_kwargs", Dict())

    LLR_min = pop!(params, "LLR_min", 0)

    props = TOML.parsefile(dirname(params["filename"]) * "/../observed_properties.toml")

    params = LilGuys.dict_to_tuple(params)
    filt_params = GaiaFilterParams(props; params...)

    stars = read_gaia_stars(filt_params)
    members = select_members(stars, filt_params)

    LilGuys.write_fits(args["output"], members)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

