#!/bin/env julia
import Random

using ArgParse
using LilGuys
using OrderedCollections



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
            default = "sculptor"
        "--time-max"
            help = "maximum time to integrate for in code units"
            default = -2120.0
            arg_type = Float64
        "--num-timesteps"
            help = "number of timesteps to record for peri / apo"
            default = 2048
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

include("orbit_utils.jl")

function main(args)
    if isfile(args["output"])
        rm(args["output"])
    end

    units = get_units(args["input"])
    pot = get_potential(args["input"])
    obs_props = get_obs_props(args["input"], args["galaxy"])

    labels, coords = get_coords(args["input"])

    orbits = LilGuys.agama_orbit(pot, coords; timerange=(0, args["time-max"]), N=args["num-timesteps"], agama_units=units)
    df_props = orbital_properties(pot, orbits, agama_units=units)

    df_icrs = LilGuys.to_frame(coords)

    gcs = LilGuys.transform.(Galactocentric, coords)
    df_galcen = LilGuys.to_frame(gcs)

    df_all = hcat(df_props, df_icrs, df_galcen)

    df_all[!, "orbit_name"] = labels

    write_fits(args["output"], df_all)

    for i in eachindex(labels)
        orbit = orbits[i]

        @assert issorted(orbit.times, rev=true)
        _, _, apos, idx_apos, _ = LilGuys.all_peris_apos(orbit)
        if length(idx_apos) > 0
            idx_i = idx_apos[end]
        else
            idx_i = length(orbit.times)
        end

        orbit_out = LilGuys.reverse(orbit[1:idx_i])

        write(joinpath(args["input"], "orbit_$(labels[i]).csv"), orbit_out)

        # save TOML
        orbit_props = OrderedDict(k => df_all[i, k] for k in names(df_all))
        orbit_props["x_i"] = orbit_out.positions[1, 1]
        orbit_props["y_i"] = orbit_out.positions[2, 1]
        orbit_props["z_i"] = orbit_out.positions[3, 1]
        orbit_props["v_x_i"] = orbit_out.velocities[1, 1]
        orbit_props["v_y_i"] = orbit_out.velocities[2, 1]
        orbit_props["v_z_i"] = orbit_out.velocities[3, 1]
        orbit_props["idx_i"] = length(orbit) - idx_i + 1
        orbit_props["t_i"] = orbit_out.times[1]
        for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
            orbit_props["$(key)_err"] = LilGuys.get_uncertainty(obs_props, key)
        end

        open(joinpath(args["input"], "orbit_$(labels[i]).toml"), "w") do f
            TOML.print(f, orbit_props)
        end
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    main(args)
end
