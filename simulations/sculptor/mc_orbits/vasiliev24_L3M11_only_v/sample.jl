import LilGuys as lguys


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs = lguys.coord_from_file(obs_props_filename)
obs = lguys.coord_err_from_file(obs_props_filename)
v_err = 6 # km/s
x_err = 1e-3 # kpc

function sample(N = 10_000)
    gc = lguys.transform(lguys.Galactocentric, obs)

    pos = hcat([lguys.position_of(gc) for _ in 1:N]...)
    vel = hcat([lguys.velocity_of(gc) for _ in 1:N]...)

    vel .+= randn(size(vel)) .* v_err
    pos .+= randn(size(pos)) .* x_err

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
    lguys.save("initial.hdf5", snap)
    println("Done")
end
