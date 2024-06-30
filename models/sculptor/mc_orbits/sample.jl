import LilGuys as lguys


const obs = lguys.ICRS(
    ra = 15.03917,
    dec = -33.70917,
    distance = 86,
    pm_ra = 0.099,
    pm_dec = -0.160,
    radial_velocity = 111.4,
)

const err = lguys.ICRS(
    ra = 0,
    dec = 0,
    distance = 3,
    pm_ra = 0.02,
    pm_dec = 0.02,
    radial_velocity = 1,
)


function sample(N = 10000)
    mc_obs = lguys.rand_coords(obs, err, N)
    pushfirst!(mc_obs, obs) # add the true observation to the sample

    mc_phase = lguys.transform.(lguys.Galactocentric, mc_obs)
    pos = hcat([p.position for p in mc_phase]...)
    vel = hcat([p.velocity for p in mc_phase]...)

    pos ./= lguys.R0
    vel ./= lguys.V0
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))
    return snap
end


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    lguys.save("initial.hdf5", snap)
end
