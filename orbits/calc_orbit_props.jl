#!/bin/env julia
import Random

using ArgParse
using LilGuys



function get_args()
    s = ArgParseSettings(
        description="""Calculates the orbits
        returing essential information.
""",
    )

    @add_arg_table s begin
        "input"
            help="input directory, should contain agama_potential.ini"
        "output"
            help="output file. defaults to orbital_properties.fits"
        "-N", "--num-orbits"
            help = "number of orbits to compute"
            default = 10_000
            arg_type=Int
        "--seed"
            help = "Random number seed. If set to -1, random seed"
            default = -1
            arg_type=Int
        "--galaxy"
            help = "name of default galaxy properties"
        "--time-max"
            help = "maximum time to integrate for in code units"
            default = -10 / T2GYR
            arg_type = Float64
        "--num-timesteps"
            help = "number of timesteps to record for peri / apo"
            default = 10_000
            arg_type = Int
    end

    args = parse_args(s)


    if isnothing(args["output"])
        args["output"] = joinpath(dirname(args["input"] * "/"), "orbital_properties.fits")
    end

    if isnothing(args["galaxy"])
        args["galaxy"] = split(args["input"], "/")[1]
    end

    return args
end


include("orbit_utils.jl")

function main(args)
    if isfile(args["output"])
        rm(args["output"])
    end

    if args["seed"] != -1
        Random.seed!(args["seed"])
    end

    units = get_units(args["input"])
    pot = get_potential(args["input"])
    obs_props = get_obs_props(args["input"], args["galaxy"])
    coords = LilGuys.rand_coords(obs_props, args["num-orbits"])

    orbits = LilGuys.agama_orbit(pot, coords; timerange=(0, args["time-max"]), N=args["num-timesteps"], agama_units=units)
    df_props = orbital_properties(pot, orbits, agama_units=units)

    df_icrs = LilGuys.to_frame(coords)

    gcs = LilGuys.transform.(Galactocentric, coords)
    df_galcen = LilGuys.to_frame(gcs)

    df_all = hcat(df_props, df_icrs, df_galcen)

    @info "writing properties"
    write_fits(args["output"], df_all)
    @info "writing orbits"
    write_orbits(args["input"], orbits)
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    main(args)
end
