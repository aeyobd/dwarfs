import LilGuys as lguys
using Distributions

import TOML

using DataFrames


function sample()
    df = TOML.parsefile("initial_conditions.toml")

    obs_special = lguys.ICRS[]
    mc_phase = lguys.Galactocentric[]

    for orbit in df["orbits"]
        println(orbit)
        icrs = icrs_from_orbit(orbit)
        frame = frame_from_orbit(orbit)

        gc = lguys.transform(lguys.Galactocentric, icrs, frame=frame)
        push!(mc_phase, gc)
    end

    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0.
    snap = lguys.Snapshot(pos, vel, m)

    return snap#, frames_df_special
end


function icrs_from_orbit(orbit)
    icrs = lguys.ICRS(
        ra = orbit["ra"],
        dec = orbit["dec"],
        distance = orbit["distance"],
        pmra = orbit["pmra"],
        pmdec = orbit["pmdec"],
        radial_velocity = orbit["radial_velocity"],
    )

    return icrs
end

function frame_from_orbit(orbit)
    if "frame" in keys(orbit)
        frame = lguys.Galactocentric(
            galcen_distance = orbit["frame"]["galcen_distance"],
            galcen_v_sun = orbit["frame"]["galcen_v_sun"],
            z_sun = orbit["frame"]["z_sun"],
            roll = orbit["frame"]["roll"],
            galcen_coord = lguys.ICRS(
                ra = orbit["frame"]["galcen_coord"]["ra"],
                dec = orbit["frame"]["galcen_coord"]["dec"],
                distance = orbit["frame"]["galcen_coord"]["distance"],
                pmra = orbit["frame"]["galcen_coord"]["pmra"],
                pmdec = orbit["frame"]["galcen_coord"]["pmdec"],
                radial_velocity = orbit["frame"]["galcen_coord"]["radial_velocity"],
            )
        )
    else
        frame = lguys.default_gc_frame
    end

    return frame
end


function (@main)(ARGS)
    snap = sample()
    # lguys.write_fits("gc_frames.fits", frames_df)
    lguys.save("initial.hdf5", snap)
end
