import LilGuys as lguys
using CSV, DataFrames


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/pace_complete_nolmc.csv"

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

    Mstar = [10 .^ o.mass_stellar ./ lguys.M2MSUN for o in eachrow(obs_props)]
    vcircmax = lguys.vel_from_M_s_fattahi.(Mstar)
    rcircmax = lguys.Ludlow.solve_rmax.(vcircmax)
    halos = [lguys.TruncNFW(v_circ_max=v, r_circ_max=r, trunc=10) for (v, r) in zip(vcircmax, rcircmax)]
    Ms = [lguys.calc_M(halo, 10) for halo in halos]
    println(Ms)

    snap = lguys.Snapshot(pos, vel, Ms)

    return snap
end


function (@main)(ARGS)
    snap = sample()
    lguys.save("initial.hdf5", snap)
end
