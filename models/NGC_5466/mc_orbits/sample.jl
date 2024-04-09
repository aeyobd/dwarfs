import LilGuys as lguys


const obs = lguys.Observation(
    ra = 211.3637,
    dec = 28.5344,
    distance = 16,
    pm_ra = -5.41,
    pm_dec = -0.79,
    radial_velocity = 106.93,
)

const err = lguys.Observation(
    ra = 0,
    dec = 0,
    distance = 0.4,
    pm_ra = 0.01,
    pm_dec = 0.01,
    radial_velocity = 0.18,
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
