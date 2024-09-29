import LilGuys as lguys


function main()
    r0 = 8
    v0 = 0.8
    r = [r0, 0, 0]
    r4 = [0, 0, 1]


    v1 = [0, v0, 0.05v0]
    v2 = [0, 1.2v0, 0.05v0]
    v3 = [0, 0.8v0, 0.05v0]
    v4 = [0, 0, 0]

    pos = hcat(r, -r, r, r4)
    vel = hcat(v1, v2, v3, v4)
    m = 0

    N = size(pos, 1)

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
