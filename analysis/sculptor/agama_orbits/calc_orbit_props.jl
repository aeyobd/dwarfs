#!/bin/env julia

import Random
import TOML
import LinearAlgebra

using DataFrames
using ArgParse

import Agama
using LilGuys
using PyFITS


Base.Broadcast.broadcastable(p::Agama.Potential) = Ref(p)
Base.Broadcast.broadcastable(p::Agama.AgamaUnits) = Ref(p)



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
            default = "sculptor"
        "--time-max"
            help = "maximum time to integrate for in code units"
            default = -10 / T2GYR
        "--num-timesteps"
            help = "number of timesteps to record for peri / apo"
            default = 10_000
    end

    args = parse_args(s)


    if isnothing(args["output"])
        args["output"] = joinpath(dirname(args["input"]), "orbital_properties.fits")
    end

    return args
end


function get_obs_props(input, galaxy)
    filename = joinpath(input, "observed_properties.toml")
    if isfile(filename)
        return TOML.parsefile(filename)
    else
        filename = joinpath(ENV["DWARFS_ROOT"], "observations/$(galaxy)/observed_properties.toml")
        obs_props = TOML.parsefile(filename)
        return obs_props
    end
end


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

    write_fits(args["output"], df_all)
end


function get_units(input)
    filename = joinpath(input, "agama_units.toml")
    if isfile(filename)
        unit_kwargs = TOML.parsefile(filename)
        return Agama.AgamaUnits(; units_kwargs...)
    else
        return Agama.AgamaUnits()
    end
end


function get_potential(directory; kwargs...)
    Agama.Potential(file=joinpath(directory, "agama_potential.ini"); kwargs...)
end


function calc_orbits(pot, coords_i; tmax=-10/T2GYR, N=10001, kwargs...)
end


function orbital_properties(pot, orbits; agama_units)
    @info "calculating peris and apos"
    peris = minimum.(radii.(orbits))
    apos = maximum.(radii.(orbits))
    
    @info "calculating orbit timescales"
    periods = orbit_period.(orbits)
    t_last = t_last_peri.(orbits)

    @info "calculating energies"
    pos = hcat((orbit.positions[:, 1] for orbit in orbits)...)
    vel = hcat((orbit.velocities[:, 1] for orbit in orbits)...)

    L = LilGuys.angular_momenta(pos, vel)
    E = energy_initial(pot, pos, vel, agama_units=agama_units)

    @info "calculating actions"
    J = actions_initial(pot, pos, vel, agama_units=agama_units)


    dt = LilGuys.mean.(LilGuys.times.(orbits))
    properties = DataFrame(
        "pericentre" => peris,
        "apocentre" => apos,
        "time_last_peri" => t_last,
        "Jr_i" => J[1, :],
        "Jz_i" => J[2, :],
        "Jphi_i" => J[3, :],
        "Lx_i" => L[1, :],
        "Ly_i" => L[2, :],
        "Lz_i" => L[3, :],
        "E_i" => E, 
        "period" => periods,
        "dt" => dt,
   )
end


function t_last_peri(orbit::LilGuys.Orbit)
    return LilGuys.last_time_peri(LilGuys.times(orbit), radii(orbit))[1]
end


function orbit_period(orbit)
    _, idxs, _, _ = LilGuys.all_peris_apos(orbit)
    period = LilGuys.mean(diff(orbit.times[idxs]))
    return period
end

to_sym_mat(x) = [x[1] x[4] x[6] 
				x[4] x[2] x[5]
				x[6] x[5] x[3]
]

function scalar_tidal_forces(pot::Agama.Potential, positions::AbstractMatrix{<:Real}; agama_units)
	T = Agama.stress(pot, positions, agama_units)
	return  maximum.(LinearAlgebra.eigvals.(to_sym_mat.(eachcol(T))))
end


function max_tidal_force(pot::Agama.Potential, orbit::LilGuys.Orbit; agama_units)
    return maximum(scalar_tidal_forces(pot, orbit.positions, agama_units))
end


function actions_initial(pot::Agama.Potential, pos, vel; agama_units)
    af = Agama.ActionFinder(pot)
    actions = Agama.actions(af, pos, vel, agama_units)
    return actions
end


function angular_momentum_initial(pos, vel)
    Ls = LilGuys.angular_momenta(pos, vel)

    return Ls
end

function energy_initial(pot::Agama.Potential, pos, vel; agama_units)
    Φ = Agama.potential.(pot, eachcol(pos), agama_units)
    T = 1/2 * radii(vel) .^ 2
    return Φ .+ T
end

if abspath(PROGRAM_FILE) == @__FILE__
    args = get_args()
    main(args)
end
