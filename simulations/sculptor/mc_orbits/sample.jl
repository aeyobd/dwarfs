import LilGuys as lguys


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)
err = lguys.coord_err_from_file(obs_props_filename)


function sample(N = 10000)
    mc_obs = lguys.rand_coords(obs, err, N)

    obs_special = [obs]
    # INPUT: add special orbits here that we actually adopt
    d_vec_1 = [-1.8, 0.87, -0.47, 0.2]
    d_vec_2 = [1.69, -0.87, 0.75, -0.05]

    for d_vec in [d_vec_1, d_vec_2]
        obs_1 = lguys.ICRS(ra = obs.ra, dec=obs.dec,
            distance = obs.distance + d_vec[1] * err.distance,
            pmra = obs.pmra + d_vec[2] * err.pmra, 
            pmdec = obs.pmdec + d_vec[3] * err.pmdec,
            radial_velocity = obs.radial_velocity + d_vec[4] * err.radial_velocity
        )

        push!(obs_special, obs_1)
    end

    mc_obs = [obs_special; mc_obs]

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


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    lguys.save("initial.hdf5", snap)
end
