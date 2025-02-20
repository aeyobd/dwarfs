import LilGuys as lguys
import TOML


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml"
obs_props = TOML.parsefile(obs_props_filename)

function sample(N = 100_000)
    mc_obs = lguys.rand_coords(obs_props, N)

    mc_phase = lguys.transform.(lguys.Galactocentric, mc_obs)
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))
    return snap
end


function (@main)(ARGS)
    println("Sampling initial conditions")
    snap = sample()
    lguys.save("initial.hdf5", snap)
    println("Done")
end
