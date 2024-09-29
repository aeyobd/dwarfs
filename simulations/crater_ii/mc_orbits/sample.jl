import LilGuys as lguys


const obs = lguys.Observation(
    ra = 177.3,
    dec = -18.4,
    distance = 117.5,
    pm_ra = -0.246,
    pm_dec = -0.227,
    radial_velocity = 87.5,
)

const obs1 = lguys.Observation(
    ra = 177.3,
    dec = -18.4,
    distance = 116,
    pm_ra = -0.169,
    pm_dec = -0.267,
    radial_velocity = 87.8,
)

const obs2 = lguys.Observation(
    ra = 177.3,
    dec = -18.4,
    distance = 117,
    pm_ra = -0.102,
    pm_dec = -0.225,
    radial_velocity = 87.2,
)


const obs3 = lguys.Observation(
    ra = 177.3,
    dec = -18.4,
    distance = 118,
    pm_ra = -0.07,
    pm_dec = -0.11,
    radial_velocity = 87.5,
)


const err = lguys.Observation(
    ra = 0,
    dec = 0,
    distance = 1.1,
    pm_ra = 0.052,
    pm_dec = 0.026,
    radial_velocity = 0.4,
)


function sample(N = 10000)
    mc_obs = lguys.rand_coords(obs, err, N)
    pushfirst!(mc_obs, obs3)
    pushfirst!(mc_obs, obs2)
    pushfirst!(mc_obs, obs1)

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
