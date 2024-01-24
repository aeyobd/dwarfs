import LilGuys as lguys


function main()
    r0 = 20
    v0 = 0.15
    r = [r0, 0, 0]

    v1 = [0, v0, 0]
    v2 = 0.5*v1
    v3 = 1.2*v1

    pos = hcat(r, -r, r)
    vel = hcat(v1, v2, v3)
    m = 0

    N = size(pos, 1)

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=zeros(N))
    lguys.save("nfw_test.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
