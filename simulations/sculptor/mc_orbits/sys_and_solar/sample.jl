import LilGuys as lguys
using Distributions

import TOML

using DataFrames

obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)
err = lguys.coord_err_from_file(obs_props_filename)
err = lguys.ICRS(
    ra=err.ra,
    dec=err.dec,
    distance=err.distance,
    pmra=0.017,
    pmdec=0.017,
    radial_velocity=err.radial_velocity
   )


mean_frame = lguys.default_gc_frame

frame_distributions = Dict{Symbol, Distribution}(
    :d => Normal(mean_frame.d, 0.033),
    :ra => Normal(mean_frame.ra, 0.0001),
    :dec => Normal(mean_frame.dec, 0.00001),
    :z_sun => Normal(mean_frame.z_sun, 0.0003),
    :v_sun_x => Normal(mean_frame.v_sun[1], 3.0),
    :v_sun_y => Normal(mean_frame.v_sun[2], 1.4),
    :v_sun_z => Normal(mean_frame.v_sun[3], 0.08),
   )

function sample(N = 100_000)

    obs_special = [obs,
        lguys.ICRS(
            ra=obs.ra,
            dec=obs.dec,
            distance = 82.458,
            pmra = 0.139,
            pmdec = -0.194,
            radial_velocity = 111.448,
           )]

    frames_special = [mean_frame, 
        lguys.GalactocentricFrame(
            d = mean_frame.d,
            ra = mean_frame.ra,
            dec = mean_frame.dec,
            z_sun = mean_frame.z_sun,
            v_sun = [12.662, 244.826, 7.798]
           )
        ]

    mc_obs = [lguys.rand_coord(obs, err) for _ in 1:N]
    frames = Vector{lguys.GalactocentricFrame}(undef, N)
    frames_df = DataFrame(Dict(
        :d => Vector{Float64}(undef, N),
        :ra => Vector{Float64}(undef, N),
        :dec => Vector{Float64}(undef, N),
        :z_sun => Vector{Float64}(undef, N),
        :v_sun_x => Vector{Float64}(undef, N),
        :v_sun_y => Vector{Float64}(undef, N),
        :v_sun_z => Vector{Float64}(undef, N)
       ))


    for i in 1:N
        for key in keys(frame_distributions)
            frames_df[i, key] = rand(frame_distributions[key])
        end

        frames[i] = lguys.GalactocentricFrame(
            d = frames_df[i, :d],
            ra = frames_df[i, :ra],
            dec = frames_df[i, :dec],
            z_sun = frames_df[i, :z_sun],
            v_sun = [
                frames_df[i, :v_sun_x],
                frames_df[i, :v_sun_y],
                frames_df[i, :v_sun_z]
            ]
        )
    end

    # add in special orbits
    mc_obs = vcat(obs_special, mc_obs)

    frames = vcat(frames_special, frames)

    frames_df_special = DataFrame(Dict(
        :d => [frame.d for frame in frames_special],
        :ra => [frame.ra for frame in frames_special],
        :dec => [frame.dec for frame in frames_special],
        :z_sun => [frame.z_sun for frame in frames_special],
        :v_sun_x => [frame.v_sun[1] for frame in frames_special],
        :v_sun_y => [frame.v_sun[2] for frame in frames_special],
        :v_sun_z => [frame.v_sun[3] for frame in frames_special]
    ))

    frames_df = vcat(frames_df_special, frames_df)

    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o, frame=frame) for (o, frame) in zip(mc_obs, frames)]
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))

    return snap, frames_df
end


function (@main)(ARGS)
    snap, frames_df = sample()
    lguys.write_fits("gc_frames.fits", frames_df)
    lguys.save("initial.hdf5", snap)
end
