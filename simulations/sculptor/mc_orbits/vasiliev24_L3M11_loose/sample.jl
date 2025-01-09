import LilGuys as lguys


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)
err = lguys.coord_err_from_file(obs_props_filename)
err_scale = 10
err = lguys.ICRS(
    ra = err_scale*err.ra,
    dec = err_scale*err.dec,
    distance = err_scale*err.distance,
    pmra = err_scale*err.pmra,
    pmdec = err_scale*err.pmdec,
    radial_velocity = err_scale*err.radial_velocity,
)

function sample(N = 10_000)
    mc_obs = lguys.rand_coords(obs, err, N)

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
