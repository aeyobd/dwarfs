import LilGuys as lguys
using CSV, DataFrames
import TOML


function sample()
    filename = joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv")

    obs_props = CSV.read(filename, DataFrame)
    obs = [lguys.ICRS(ra=o.ra, dec=o.dec, distance=o.distance, pmra=o.pmra, pmdec=o.pmdec, radial_velocity=o.radial_velocity) for o in eachrow(obs_props)]
    index = obs_props.index

    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o) for o in obs]
    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0.0
    snap = lguys.Snapshot(pos, vel, m, index=index)

    return snap
end


function (@main)(ARGS)
    snap = sample()
    lguys.write("initial.hdf5", snap)
end
