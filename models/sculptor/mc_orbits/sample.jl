import LilGuys as lguys


const obs = lguys.ICRS(
    ra = 15.0183,
    dec = -33.7186,
    distance = 83.2,
    pmra = 0.099,
    pmdec = -0.160,
    radial_velocity = 111.03,
)

const err = lguys.ICRS(
    ra = 0.02, #0.0012,
    dec = 0.01, #0.00072,
    distance = 2,
    pmra = 0.02,
    pmdec = 0.02,
    radial_velocity = 0.3, #0.2,
)


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
