import LilGuys as lguys
import TOML


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"
obs_props = TOML.parsefile(obs_props_filename)

obs = lguys.coord_from_file(obs_props_filename)
err = lguys.ICRS(;Dict(sym => 2*lguys.get_uncertainty(obs_props, string(sym)) for sym in [:ra, :dec, :distance, :pmra, :pmdec, :radial_velocity])...)


function sample(N = 100_000)
    mc_obs = lguys.rand_coords(obs, err, N)

    mc_phase = lguys.transform.(lguys.Galactocentric, mc_obs)
    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))
    return snap
end


function (@main)(ARGS)
    println("Sampling initial conditions")
    snap = sample()
    lguys.write("initial.hdf5", snap)
    println("Done")
end
