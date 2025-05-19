import LilGuys as lguys
using Distributions

import TOML

using DataFrames

obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)


function sample(N = 100_000)
    obs_special = [obs,
        lguys.ICRS(
            ra=obs.ra,
            dec=obs.dec,
            distance = 86.21,
            pmra = 0.14,
            pmdec = -0.15,
            radial_velocity = 111.40,
           ),
        lguys.ICRS(
            ra=obs.ra,
            dec=obs.dec,
            distance = 84.98,
            pmra = 0.06,
            pmdec = -0.13,
            radial_velocity = 111.42,
           ),
       ]


    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o) for o in obs_special]
    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0.
    snap = lguys.Snapshot(pos, vel, m)

    return snap
end


function (@main)(ARGS)
    snap = sample()
    lguys.write("initial.hdf5", snap)
end
