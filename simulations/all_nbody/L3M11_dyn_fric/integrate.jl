import Agama
import LilGuys
import Random
import LilGuys: radii

Random.seed!(127)

module NBody 
	include("../nbody_utils.jl")
end




function random_orbit(pot, obs_props; units, timestep=0.05, f_σ, f_ρ)
    ic = NBody.sample_initial_conditions(obs_props)

    traj_evolving = NBody.integrate_particles(
        ic...,
        NBody.force_potential_nbody_friction(pot, units, ic[3], f_σ=f_σ, f_ρ=f_ρ),
        timestep=timestep,
    )
end


"""
    write_orbits(output, orbits; N_max)

Write the first `N_max` orbits to a file "orbits.hdf5" in `output`.
"""
function write_orbits(filename, orbits; N_max=1000)

    N_max = min(length(orbits), N_max)
    structs = [(string(i) => orbit) for (i, orbit) in enumerate(orbits[1:N_max])]

    LilGuys.write_structs_to_hdf5(filename, structs)
end





function (@main)(ARGS)
    pot = NBody.get_potential_c("vasiliev24/L3M11/potential")
    units = Agama.VASILIEV_UNITS
    obs_props = NBody.get_obs_props()
    Norbits = 10

    pot_sigma = NBody.get_potential("vasiliev24/L3M11/potential_mw_init")
    σv = NBody.get_sigma_v(pot_sigma, units)
    f_σ(x) = σv(radii(x))
    ρ = NBody.get_rho(pot_sigma, units)
    f_ρ(x) = ρ(radii(x))


    if !isdir("out")
        mkdir("out")
    end

    for i in 1:Norbits
        orbits = random_orbit(pot, obs_props, units=units, f_σ=f_σ, f_ρ=f_ρ)
        write_orbits("out/orbits_$i.hdf5", orbits)
    end
end
