import LilGuys as lguys
using CSV, DataFrames


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/all_galaxies.csv"

obs_props = CSV.read(obs_props_filename, DataFrame)
filt = .!isnan.(obs_props.radial_velocity)
obs_props = obs_props[filt, :]


function sample()
    obs= [lguys.ICRS(ra=o.ra, dec=o.dec, distance=o.distance, pmra=o.pmra, pmdec=o.pmdec, radial_velocity=o.radial_velocity) for o in eachrow(obs_props)]


    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o) for o in obs]
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, size(pos, 2)))

    return snap
end


function (@main)(ARGS)
    snap = sample()
    lguys.save("initial.hdf5", snap)
end
