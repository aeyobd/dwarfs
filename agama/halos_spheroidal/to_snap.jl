using LilGuys
using CSV, DataFrames



function (@main)(ARGS)
    filename_in = ARGS[1]
    outname = splitext(filename_in)[1] * ".hdf5"

    df = CSV.read(filename_in, DataFrame, comment="#", header=["x", "y", "z", "v_x", "v_y", "v_z", "mass"])
    pos = [df.x df.y df.z]'
    vel = [df.v_x df.v_y df.v_z]'
    mass = df.mass

    snap = LilGuys.Snapshot(pos, vel, mass)
    LilGuys.write(outname, snap)
end
