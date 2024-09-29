import LilGuys as lguys


obs_prop_filename = "../../../sculptor_obs_properties.toml"
obs = lguys.coord_from_file(obs_props_filename)
err = lguys.coord_err_from_file(obs_props_filename)


function sample(N = 10000)
    mc_obs = lguys.rand_coords(obs, err, N)
    pushfirst!(mc_obs, obs) # add the true observation to the sample

    mc_phase = lguys.transform.(lguys.Galactocentric, mc_obs)
    pos = hcat([p.position for p in mc_phase]...)
    vel = hcat([p.velocity for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))
    return snap
end


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    lguys.save("initial.hdf5", snap)
end
