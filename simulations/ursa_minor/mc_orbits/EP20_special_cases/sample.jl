import LilGuys as lguys
using Distributions

import TOML

using DataFrames

obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/ursa_minor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)


function sample()
    mean_frame = lguys.default_gc_frame

    obs_special = [obs,
        lguys.ICRS(
            ra=obs.ra,
            dec=obs.dec,
            distance = 65.917,
            pmra = -0.160,
            pmdec = 0.051,
            radial_velocity = -245.258,
           ),
        lguys.ICRS(
            ra = obs.ra,
            dec = obs.dec,
            distance = 75.172,
            pmra = -0.093,
            pmdec = 0.107,
            radial_velocity = -245.159,
           )
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

    m = 0
    snap = lguys.Snapshot(pos, vel, m)

    return snap#, frames_df_special
end


function (@main)(ARGS)
    snap = sample()
    # lguys.write_fits("gc_frames.fits", frames_df)
    lguys.save("initial.hdf5", snap)
end
