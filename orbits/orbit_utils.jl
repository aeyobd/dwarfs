import TOML
import LinearAlgebra

import Agama
using LilGuys

using DataFrames, CSV
using PyFITS


Base.Broadcast.broadcastable(p::Agama.Potential) = Ref(p)
Base.Broadcast.broadcastable(p::Agama.AgamaUnits) = Ref(p)


"""
    get_obs_props(orbit_dir, galaxy)

Retrieves the observed properties, either from the directory `orbit_dir/observed_properties.toml` or 
the defaults for that galaxy in the `/observations` directory.
"""
function get_obs_props(orbit_dir::String, galaxy::String)
    filename = joinpath(orbit_dir, "observed_properties.toml")
    if isfile(filename)
        @info "loading $filename"
        return TOML.parsefile(filename)
    else
        filename = joinpath(ENV["DWARFS_ROOT"], "observations/$(galaxy)/observed_properties.toml")
        obs_props = TOML.parsefile(filename)
        return obs_props
    end
end

"""
    get_units(orbit_dir)

Retrieve the units specified in the orbit model directory, or returns our default
scalefree units (as an `Agama.AgamaUnits()` object).
"""
function get_units(orbit_dir::String)
    filename = joinpath(orbit_dir, "agama_units.toml")
    if isfile(filename)
        unit_kwargs = TOML.parsefile(filename) |> LilGuys.dict_to_tuple
        return Agama.AgamaUnits(; unit_kwargs...)
    else
        return Agama.AgamaUnits()
    end
end


"""
    get_potential(orbit_dir; kwargs...)

Retrieve the potential from the specified directory. kwargs passed to Agama.
"""
function get_potential(orbit_dir::String; kwargs...)
    Agama.Potential(file=joinpath(orbit_dir, "agama_potential.ini"); kwargs...)
end


"""
    get_initial_apocentre(orbit)

Find the apocentre occuring earliest in time for the given orbit, used
for our initial conditions.
"""
function get_initial_apocentre(orbit)
    @assert issorted(orbit.times, rev=true)

    _, _, apos, idx_apos, _ = LilGuys.all_peris_apos(orbit)
    if length(idx_apos) > 0
        idx_i = idx_apos[end]
    else
        idx_i = length(orbit.times)
    end

    return idx_i
end
     

function orbital_properties(pot, orbits; agama_units, actions=false)
    @assert all([issorted(orbit.times, rev=true) for orbit in orbits])

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


    dt = [LilGuys.mean(diff(LilGuys.times(orbit))) for orbit in orbits]
    properties = DataFrame(
        "pericentre" => peris,
        "apocentre" => apos,
        "n_peris" => n_peris,
        "period" => periods,
        "period_apo" => period_apos,
        "time_last_peri" => t_last,
        "peri_apo_max_err" => errs,
        "orbital_timestep" => dt,

        "x_i" => [orbit.positions[1, end] for orbit in orbits],
        "y_i" => [orbit.positions[2, end] for orbit in orbits],
        "z_i" => [orbit.positions[3, end] for orbit in orbits],
        "v_x_i" => [orbit.velocities[1, end] for orbit in orbits],
        "v_y_i" => [orbit.velocities[2, end] for orbit in orbits],
        "v_z_i" => [orbit.velocities[3, end] for orbit in orbits],
        "t_i" => [orbit.times[end] for orbit in orbits],

        "x" => [orbit.positions[1, 1] for orbit in orbits],
        "y" => [orbit.positions[2, 1] for orbit in orbits],
        "z" => [orbit.positions[3, 1] for orbit in orbits],
        "v_x" => [orbit.velocities[1, 1] for orbit in orbits],
        "v_y" => [orbit.velocities[2, 1] for orbit in orbits],
        "v_z" => [orbit.velocities[3, 1] for orbit in orbits],
        "t_f" => [orbit.times[1] for orbit in orbits],

        "Lx_f" => L[1, :],
        "Ly_f" => L[2, :],
        "Lz_f" => L[3, :],
        "E_f" => E, 
   )


    if actions
        properties[!, "Jr_i"] = J[1, :]
        properties[!, "Jz_i"] = J[2, :]
        properties[!, "Jphi_i"] = J[3, :]
    end

    properties
end


"""
    t_last_peri(orbit)

Retrieve the time of the last pericentre for the given orbit.
"""
function t_last_peri(orbit::LilGuys.Orbit)
    return LilGuys.last_time_peri(LilGuys.times(orbit), radii(orbit))[1]
end


"""
    orbit_period(orbit)

Compute the orbital period and related properties for an orbit. 
Returns a tuple of the period between pericentres, the estimated uncertainty
of this period, the period between apocentres, and the number of pericentres.
"""
function orbit_period(orbit)
    peris, idxs, apos, idx_apos, err = LilGuys.all_peris_apos(orbit)
    period = abs.(LilGuys.mean(diff(orbit.times[idxs])))
    if length(idx_apos) > 1
        period_apo = abs(diff(orbit.times[idx_apos])[1]) # since most recent is lower index
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


"""
    scalar_tidal_forces(pot, orbit; agama_units)

Calculate the scalar tidal force (max eigenvalue of the tidal strain matrix) for 
each position along the orbit. `orbit` may be given as a matrix of positions or a 
`LilGuys.Orbit` object.
"""
function scalar_tidal_forces(pot::Agama.Potential, orbit::LilGuys.Orbit; agama_units)
    return scalar_tidal_forces(pot, orbit.positions; agama_units=agama_units, t=orbit.times)
end


"""
    max_tidal_force(pot, orbit; agama_units)

Given a potential and a orbit (or a matrix of positions), 
calculate the maximum tidal force along the orbit. 
See `scalar_tidal_forces` as well.
"""
function max_tidal_force(pot::Agama.Potential, orbit; agama_units)
    return maximum(scalar_tidal_forces(pot, orbit; agama_units=agama_units))
end


"""
    actions_initial(pot, pos, vel; agama_units)

Given the potential and initial positions and velocities, returns the initial actions.
"""
function actions_initial(pot::Agama.Potential, pos, vel; agama_units)
    af = Agama.ActionFinder(pot)
    actions = Agama.actions(af, pos, vel, agama_units)
    return actions
end


"""
    energy_initial(pot::Agama.Potential, pos::Matrix, vel::Matrix; agama_units)

Return the initial energy for particles at a given position and velocity
"""
function energy_initial(pot::Agama.Potential, pos::AbstractMatrix{<:Real}, vel::AbstractMatrix{<:Real}; agama_units)
    Φ = Agama.potential.(pot, eachcol(pos), agama_units)
    T = 1/2 * radii(vel) .^ 2
    return Φ .+ T
end


"""
    write_orbits(output_dir, orbits; N_max)

Write the first `N_max` orbits to a file "orbits.hdf5" in `output_dir`.
"""
function write_orbits(output_dir, orbits; N_max=1000)
    filename = joinpath(output_dir, "orbits.hdf5")

    N_max = min(length(orbits), N_max)
    structs = [(string(i) => orbit) for (i, orbit) in enumerate(orbits[1:N_max])]

    LilGuys.write_structs_to_hdf5(filename, structs)
end




"""
    get_lmc_orbit(orbit_dir)

Retrieve the orbit of the LMC as a LilGuys.Orbit object assuming a vasiliev+21 like potential.
"""
function get_lmc_orbit(orbit_dir::String)
    lmc_file = joinpath(orbit_dir, "trajlmc.txt")
    df_lmc = lmc_traj = CSV.read(lmc_file, DataFrame, delim=" ", header = [:time, :x, :y, :z, :v_x, :v_y, :v_z], ignorerepeated=true, ntasks=1)

    pos = hcat(df_lmc.x, df_lmc.y, df_lmc.z)'
    vel = hcat(df_lmc.v_x, df_lmc.v_y, df_lmc.v_z)'

    # convert to code units
    t = df_lmc.time .* Agama.time_scale(Agama.VASILIEV_UNITS) 
    vel .*= Agama.velocity_scale(Agama.VASILIEV_UNITS) 

    orbit_lmc = Orbit(times=t, positions=pos, velocities=vel)
end
