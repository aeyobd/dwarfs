using FITSIO
import LilGuys as lguys

function fits_to_snap(fitsfile::String, snapfile::String = "")
    pos = []
    vel = []
    mass = []

    println("Reading data from $fitsfile")
    FITS(fitsfile, "r") do f
        data = f[2]

        pos = [read(data, "x") read(data, "y") read(data, "z")]'
        vel = [read(data, "vx") read(data, "vy") read(data, "vz")]'
        mass = read(data, "mass")
        
    end

    println("Creating snapshot")
    snap = lguys.Snapshot(pos, vel, mass)

    if snapfile == ""
        snapfile = splitext(fitsfile)[1] * ".hdf5"
    end

    println("Writing snapshot")
    lguys.save(snapfile, snap)
end



function read_args()
    if length(ARGS) < 1
        println("Usage: fits_to_snap.jl <fitsfile> [snapfile]")
        exit(1)
    end

    fitsfile = ARGS[1]
    snapfile = ""
    if length(ARGS) > 1
        snapfile = ARGS[2]
    end

    return fitsfile, snapfile
end


if abspath(PROGRAM_FILE) == @__FILE__
    fitsfile, snapfile = read_args()
    fits_to_snap(fitsfile, snapfile)
end
