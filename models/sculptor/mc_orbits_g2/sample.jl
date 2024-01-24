using LilGuys


const obs = LilGuys.Observation(
    ra = 15.03917,
    dec = -33.70917,
    distance = 86,
    pm_ra = 0.099,
    pm_dec = -0.160,
    radial_velocity = 111.4,
)

const err = LilGuys.Observation(
    ra = 0.00001,
    dec = 0.00001,
    distance = 5,
    pm_ra = 0.02,
    pm_dec = 0.02,
    radial_velocity = 3
)


function sample(N = 1000)
    mc_obs = LilGuys.rand_coords(obs, err, N)
    pushfirst!(mc_obs, obs) # add the true observation to the sample
    mc_phase = [LilGuys.to_galcen(o) for o in mc_obs]
    pos = reduce(hcat, [[p.position...] for p in mc_phase])
    vel = reduce(hcat, [[p.velocity...] for p in mc_phase])

    vel *= 1 # reverse velocities to go backwards in time
    m = 0
    header = LilGuys.make_gadget2_header(N, m)
    snap = Snapshot(positions=pos, velocities=vel, header=header, masses=LilGuys.ConstVector(m, N), index=collect(1:N))
    return snap
end


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    save("positions.hdf5", snap)
end
