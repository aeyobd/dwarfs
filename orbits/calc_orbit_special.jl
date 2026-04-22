#!/bin/env julia
import Random

using ArgParse
using LilGuys
using OrderedCollections


include("orbit_utils.jl")

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
        "--galaxy"
            help = "name of default galaxy properties"
        "--time-max"
            help = "maximum time to integrate for in code units"
            default = -2120.0
            arg_type = Float64
        "--num-timesteps"
            help = "number of timesteps to record for peri / apo"
            default = 10_000
            arg_type = Int
    end

    args = parse_args(s)


    if isnothing(args["output"])
        args["output"] = joinpath(dirname(args["input"]), "orbital_properties.fits")
    end

    return args
end


function get_coords(input)
    params = TOML.parsefile(joinpath(input, "initial_conditions.toml"))

    labels = String[]
    coords = ICRS[]
    for row in params["orbits"]
        push!(labels, row["name"])
        push!(coords, ICRS(row))
    end

    return labels, coords
end

function get_time_max(input, time_max_default)
    params = TOML.parsefile(joinpath(input, "initial_conditions.toml"))

    time_maxes = Float64[]
    for row in params["orbits"]
        t = get(row, "time_max", time_max_default)
        push!(time_maxes, t)
    end

    return time_maxes
end


function main(args)
    if isfile(args["output"])
        rm(args["output"])
    end

    agama_units = get_units(args["input"])
    pot = get_potential(args["input"])
    obs_props = get_obs_props(args["input"], args["galaxy"])

    labels, coords = get_coords(args["input"])
    time_maxes = get_time_max(args["input"], args["time-max"])
    timerange = [(0, time_max) for time_max in time_maxes]

    orbits = LilGuys.agama_orbit(pot, coords; timerange=timerange, 
                                 N=args["num-timesteps"], agama_units=agama_units)
    df_props = orbital_properties(pot, orbits, agama_units=agama_units)

    df_icrs = LilGuys.to_frame(coords)

    df_all = hcat(df_props, df_icrs)

    df_all[!, "orbit_name"] = labels

    write_fits(args["output"], df_all)

    write_individual_orbits(labels, coords, orbits; agama_units=agama_units, obs_props=obs_props, pot=pot)
end

function write_individual_orbits(labels, coords, orbits; agama_units, obs_props, pot)
    for i in eachindex(labels)
        idx_i = get_initial_apocentre(orbits[i])
        orbit = orbits[i][1:idx_i]

        write(joinpath(args["input"], "orbit_$(labels[i]).csv"), LilGuys.reverse(orbit))

        # recalculate on truncated orbit
        df_props = orbital_properties(pot, [orbit], agama_units=agama_units)
        orbit_props = OrderedDict{String, Any}(k => df_props[1, k] for k in names(df_props))
        orbit_props["idx_i"] = length(orbit) - idx_i + 1

        # add observation errors
        for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
            orbit_props["$(key)"] = getproperty(coords[i], Symbol(key))
            orbit_props["$(key)_err"] = LilGuys.get_uncertainty(obs_props, key)
        end
        orbit_props["orbit_name"] = labels[i]

        # save TOML
        open(joinpath(args["input"], "orbit_$(labels[i]).toml"), "w") do f
            TOML.print(f, orbit_props)
        end
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    main(args)
end
