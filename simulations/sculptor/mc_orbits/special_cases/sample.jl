import LilGuys as lguys
using Distributions

import TOML

using DataFrames

obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)


function sample(N = 100_000)
    mean_frame = lguys.default_gc_frame

    obs_special = [obs,
        lguys.ICRS(
            ra=obs.ra,
            dec=obs.dec,
            distance = 82.7,
            pmra = 0.134,
            pmdec = -0.196,
            radial_velocity = 111.4,
           ),
        lguys.ICRS(
            ra = obs.ra,
            dec = obs.dec,
            distance = 85.4,
            pmra = 0.064,
            pmdec = -0.126,
            radial_velocity = 111.4,
           )
       ]

    frames_special = [mean_frame, mean_frame, mean_frame]
    # frames_special = [mean_frame, 
    #     lguys.GalactocentricFrame(
    #         d = mean_frame.d,
    #         ra = mean_frame.ra,
    #         dec = mean_frame.dec,
    #         z_sun = mean_frame.z_sun,
    #         v_sun = [13.486, 246.457, 7.781]
    #        ),
    #     lguys.GalactocentricFrame(
    #         d = mean_frame.d,
    #         ra = mean_frame.ra,
    #         dec = mean_frame.dec,
    #         z_sun = mean_frame.z_sun,
    #         v_sun = [12.171, 244.604, 7.775]
    #        )
    #     ]

    # frames_df_special = DataFrame(Dict(
    #     :d => [frame.d for frame in frames_special],
    #     :ra => [frame.ra for frame in frames_special],
    #     :dec => [frame.dec for frame in frames_special],
    #     :z_sun => [frame.z_sun for frame in frames_special],
    #     :v_sun_x => [frame.v_sun[1] for frame in frames_special],
    #     :v_sun_y => [frame.v_sun[2] for frame in frames_special],
    #     :v_sun_z => [frame.v_sun[3] for frame in frames_special]
    # ))


    mc_obs = obs_special
    frames = frames_special
    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o, frame=frame) for (o, frame) in zip(mc_obs, frames)]
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))

    return snap#, frames_df_special
end


function (@main)(ARGS)
    snap = sample()
    # lguys.write_fits("gc_frames.fits", frames_df)
    lguys.save("initial.hdf5", snap)
end
