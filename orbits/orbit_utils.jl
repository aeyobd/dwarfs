import TOML
import LinearAlgebra

using DataFrames, CSV

import Agama
using LilGuys
using PyFITS


Base.Broadcast.broadcastable(p::Agama.Potential) = Ref(p)
Base.Broadcast.broadcastable(p::Agama.AgamaUnits) = Ref(p)


function get_obs_props(input, galaxy)
    filename = joinpath(input, "observed_properties.toml")
    if isfile(filename)
        @info "loading $filename"
        return TOML.parsefile(filename)
    else
        filename = joinpath(ENV["DWARFS_ROOT"], "observations/$(galaxy)/observed_properties.toml")
        obs_props = TOML.parsefile(filename)
        return obs_props
    end
end


function get_units(input)
    filename = joinpath(input, "agama_units.toml")
    if isfile(filename)
        unit_kwargs = TOML.parsefile(filename) |> LilGuys.dict_to_tuple
        return Agama.AgamaUnits(; unit_kwargs...)
    else
        return Agama.AgamaUnits()
    end
end


function get_potential(directory; kwargs...)
    Agama.Potential(file=joinpath(directory, "agama_potential.ini"); kwargs...)
end


function calc_orbits(pot, coords_i; tmax=-10/T2GYR, N=10001, kwargs...)
end


function orbital_properties(pot, orbits; agama_units, actions=false)
    @info "calculating peris and apos"
    peris = minimum.(radii.(orbits))
    apos = maximum.(radii.(orbits))
    
    @info "calculating orbit timescales"
    periods_errs = orbit_period.(orbits)
    periods = [p[1] for p in periods_errs]
    errs = [p[2] for p in periods_errs]
    period_apos = [p[3] for p in periods_errs]
    n_peris = [p[4] for p in periods_errs]

    max_err = maximum(errs)
    @info "max pericentre err: $max_err"

    t_last = t_last_peri.(orbits)

    @info "calculating energies"
    pos = hcat((orbit.positions[:, 1] for orbit in orbits)...)
    vel = hcat((orbit.velocities[:, 1] for orbit in orbits)...)

    L = LilGuys.angular_momenta(pos, vel)
    E = energy_initial(pot, pos, vel, agama_units=agama_units)

    if actions
        @info "calculating actions"
        J = actions_initial(pot, pos, vel, agama_units=agama_units)
    end


    dt = LilGuys.mean.(LilGuys.times.(orbits))
    properties = DataFrame(
        "pericentre" => peris,
        "apocentre" => apos,
        "n_peris" => n_peris,
        "time_last_peri" => t_last,
        "Lx_f" => L[1, :],
        "Ly_f" => L[2, :],
        "Lz_f" => L[3, :],
        "E_f" => E, 
        "period" => periods,
        "period_apo" => period_apos,
        "dt" => dt,
        "peri_apo_max_err" => errs,
        "x_i" => [orbit.positions[1, end] for orbit in orbits],
        "y_i" => [orbit.positions[2, end] for orbit in orbits],
        "z_i" => [orbit.positions[3, end] for orbit in orbits],
        "v_x_i" => [orbit.velocities[1, end] * V2KMS for orbit in orbits],
        "v_y_i" => [orbit.velocities[1, end] * V2KMS for orbit in orbits],
        "v_z_i" => [orbit.velocities[1, end] * V2KMS for orbit in orbits],
   )


    if actions
        properties[!, "Jr_i"] = J[1, :]
        properties[!, "Jz_i"] = J[2, :]
        properties[!, "Jphi_i"] = J[3, :]
    end

    properties
end


function t_last_peri(orbit::LilGuys.Orbit)
    return LilGuys.last_time_peri(LilGuys.times(orbit), radii(orbit))[1]
end


function orbit_period(orbit)
    peris, idxs, apos, idx_apos, err = LilGuys.all_peris_apos(orbit)
    period = LilGuys.mean(diff(orbit.times[idxs]))
    if length(idx_apos) > 1
        period_apo = diff(orbit.times[idx_apos])[1] # since most recent is lower index
    else
        period_apo = NaN
    end

    n_peris = length(idxs)
    return period, err, period_apo, n_peris
end


to_sym_mat(x) = [x[1] x[4] x[6] 
				x[4] x[2] x[5]
				x[6] x[5] x[3]
]

function scalar_tidal_forces(pot::Agama.Potential, positions::AbstractMatrix{<:Real}; agama_units, t=0)
	T = Agama.stress(pot, positions, agama_units, t=t)
	return  maximum.(LinearAlgebra.eigvals.(to_sym_mat.(eachcol(T))))
end


function scalar_tidal_forces(pot::Agama.Potential, orbit::LilGuys.Orbit; agama_units)
    return scalar_tidal_forces(pot, orbit.positions; agama_units=agama_units, t=orbit.times)
end


function max_tidal_force(pot::Agama.Potential, orbit::LilGuys.Orbit; agama_units)
    return maximum(scalar_tidal_forces(pot, orbit; agama_units=agama_units))
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


"""
    write_orbits(output, orbits; N_max)

Write the first `N_max` orbits to a file "orbits.hdf5" in `output`.
"""
function write_orbits(output, orbits; N_max=1000)
    filename = joinpath(output, "orbits.hdf5")

    N_max = min(length(orbits), N_max)
    structs = [(string(i) => orbit) for (i, orbit) in enumerate(orbits[1:N_max])]

    LilGuys.write_structs_to_hdf5(filename, structs)
end




function get_lmc_orbit(input)
    lmc_file = joinpath(input, "trajlmc.txt")
    df_lmc = lmc_traj = CSV.read(lmc_file, DataFrame, delim=" ", header = [:time, :x, :y, :z, :v_x, :v_y, :v_z], ignorerepeated=true, ntasks=1)

    pos = hcat(df_lmc.x, df_lmc.y, df_lmc.z)'
    vel = hcat(df_lmc.v_x, df_lmc.v_y, df_lmc.v_z)'

    # convert to code units
    t = df_lmc.time .* Agama.time_scale(Agama.VASILIEV_UNITS) 
    vel .*= Agama.velocity_scale(Agama.VASILIEV_UNITS) 

    orbit_lmc = Orbit(times=t, positions=pos, velocities=vel)
end
