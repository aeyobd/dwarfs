#!/bin/env julia

using ArgParse
using LilGuys
import PythonCall # for fits
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

    profile_kwargs = pop!(params, "profile_kwargs", Dict())

    LLR_min = pop!(params, "LLR_min", 0)

    props = TOML.parsefile(dirname(params["filename"]) * "/../observed_properties.toml")

    params = LilGuys.dict_to_tuple(params)
    filt_params = GaiaFilterParams(props; params...)

    stars = read_gaia_stars(filt_params)
    members = select_members(stars, filt_params)

    if "bins" ∈ keys(params)
        if profile_kwargs["bins"] == "equal-number"
            profile_kwargs["bins"] = LilGuys.bins_equal_number
        end
    else
        num_per_bin = ceil(Int64, LilGuys.default_n_per_bin(members.R_ell))
        num_per_bin = pop!(profile_kwargs, "num_per_bin", num_per_bin)
        bin_width = pop!(profile_kwargs, "bin_width", 0.05)
        profile_kwargs["bins"] = (x, w) -> LilGuys.bins_both(x, w, num_per_bin=num_per_bin, bin_width=bin_width)
    end


    profile_kwargs = LilGuys.dict_to_tuple(profile_kwargs)
    prof = LilGuys.StellarDensityProfile(members.R_ell; R_units="arcmin", profile_kwargs...)

    open(args["output"], "w") do f
        print(f, prof)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

