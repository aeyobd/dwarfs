import LilGuys as lguys

function make_init()
    println("loading snap")
    snap = lguys.Snapshot("initial_g4.hdf5")

    N = length(snap)
    m = 0
    snap.header = lguys.make_gadget2_header(N, m)

    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    make_init()
end
