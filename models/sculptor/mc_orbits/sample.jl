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
    ra = 0,
    dec = 0,
    distance = 3,
    pm_ra = 0.002,
    pm_dec = 0.002,
    radial_velocity = 0.37
)


function sample(N = 10_000)
    mc_obs = LilGuys.rand_coords(obs, err, N)
    mc_phase = [LilGuys.to_galcen(o) for o in mc_obs]
    pos = reduce(hcat, [[p.pos...] for p in mc_phase])
    vel = reduce(hcat, [[p.vel...] for p in mc_phase])

    vel *= -1 # reverse velocities to go backwards in time
    println(size(pos))
    println(size(vel))
    snap = Snapshot(positions=pos, velocities= vel, masses=LilGuys.ConstVector(0., N))
    return snap
end


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    save("positions.hdf5", snap)
end
