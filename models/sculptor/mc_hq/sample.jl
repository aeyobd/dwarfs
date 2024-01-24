using LilGuys


const obs = LilGuys.Observation(
    ra = 15.03917,
    dec = 33.70917,
    distance = 100, #86,
    pm_ra = 0.099,
    pm_dec = -0.160,
    radial_velocity = 111.4,
)

const err = LilGuys.Observation(
    ra = 0.1,
    dec = 0.1,
    distance = 3,
    pm_ra = 0.02,
    pm_dec = 0.02,
    radial_velocity = 10,
)


function sample(N = 1)
    mc_obs = LilGuys.rand_coords(obs, err, N)
    pushfirst!(mc_obs, obs) # add the true observation to the sample

    mc_phase = [LilGuys.to_galcen(o) for o in mc_obs]
    pos = reduce(hcat, [[p.position...] for p in mc_phase])
    vel = reduce(hcat, [[p.velocity...] for p in mc_phase])

    vel *= -1 # reverse velocities to go backwards in time
    m = 1e-20
    snap = Snapshot(pos, vel, LilGuys.ConstVector(m, N+1))
    return snap
end


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    save("positions.hdf5", snap)
end
