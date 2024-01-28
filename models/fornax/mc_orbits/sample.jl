import LilGuys as lguys


const obs = lguys.Observation(
    ra = 39.99708333,
    dec = -34.44916667,
    distance = 147,
    pm_ra = 0.374,
    pm_dec = -0.401,
    radial_velocity = 55.3,
)

const err = lguys.Observation(
    ra = 0,
    dec = 0,
    distance = 12,
    pm_ra = 0.035,
    pm_dec = 0.035,
    radial_velocity = 0.3,
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
