#!/bin/env julia

using ArgParse
using LilGuys
import TOML
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")


function get_args()

    s = ArgParseSettings(description="applies the given filters to a data file and computes the stellar profile")


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
        args["output"] = "$(filebase)_profile.toml"
    end

    return args
end

function main()
    args = get_args()

    params = read_paramfile(args["input"])
    if "profile_kwargs" ∈ keys(params)
        kwargs = pop!(params, "profile_kwargs")
    else
        println("No profile_kwargs found in input file, using default")
        kwargs = Dict()
    end

    if "bins" ∈ keys(params)
        if kwargs["bins"] == "equal-number"
            kwargs["bins"] = LilGuys.bins_equal_number
        end
    else
        kwargs["bins"] = (x, w) -> LilGuys.bins_both(x, w, num_per_bin=20, bin_width=0.05)
    end

    kwargs = LilGuys.dict_to_tuple(kwargs)

    filt_params = GaiaFilterParams(params)

    # calculates r_ell in arcminutes
    stars = read_gaia_stars(filt_params.filename, filt_params)
    members = select_members(stars, filt_params)

    prof = LilGuys.StellarProfile(members.r_ell; r_units="arcmin", kwargs...)

    open(args["output"], "w") do f
        print(f, prof)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

