import LilGuys as lguys
using CSV, DataFrames
import TOML


function sample()
    filename = joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv")

    obs_props = CSV.read(filename, DataFrame)

    obs = [lguys.ICRS(ra=o.ra, dec=o.dec, distance=o.distance, pmra=o.pmra, pmdec=o.pmdec, radial_velocity=o.radial_velocity) for o in eachrow(obs_props)]
    index = obs_props.index

    M_L_star = 2.0
    L = @. lguys.mag_to_L(obs_props.Mv)
    Mstar = @. L / lguys.M2MSUN * M_L_star

    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o) for o in obs]
    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    vcircmax = lguys.vel_from_M_s_fattahi.(Mstar)
    rcircmax = lguys.Ludlow.solve_rmax.(vcircmax)
    halos = [lguys.TruncNFW(v_circ_max=v, r_circ_max=r, trunc=10) for (v, r) in zip(vcircmax, rcircmax)]
    Ms = [lguys.M200(halo, 10) for halo in halos]


    df = copy(obs_props)
    df[!, :v_circ_max] = vcircmax
    df[!, :r_circ_max] = rcircmax
    df[!, :masses] = Ms

    CSV.write("initial_galaxies.csv", df)
    snap = lguys.Snapshot(pos, vel, Ms, index=index)

    return snap
end


function (@main)(ARGS)
    snap = sample()
    lguys.write("initial.hdf5", snap)
end
