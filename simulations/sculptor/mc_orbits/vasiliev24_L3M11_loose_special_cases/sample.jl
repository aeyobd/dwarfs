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
            ra=15.01,
            dec=-33.72,
            distance = 57.01,
            pmra = 0.43,
            pmdec = -0.008,
            radial_velocity = 111.38,
           ),
        lguys.ICRS(
            ra=15.03,
            dec=-33.72,
            distance = 141.74,
            pmra = -0.030,
            pmdec = -0.041,
            radial_velocity = 110.85,
           ),
       ]

    frames_special = [mean_frame, mean_frame, mean_frame]

    mc_obs = obs_special
    frames = frames_special
    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o, frame=frame) for (o, frame) in zip(mc_obs, frames)]
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0.
    snap = lguys.Snapshot(pos, vel, m)

    return snap#, frames_df_special
end


function (@main)(ARGS)
    snap = sample()
    # lguys.write_fits("gc_frames.fits", frames_df)
    lguys.save("initial.hdf5", snap)
end
