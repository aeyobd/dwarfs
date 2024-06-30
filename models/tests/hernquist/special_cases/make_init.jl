import LilGuys as lguys


function main()
    r0 = 2
    v0 = 0.5
    r = [r0, 0, 0]

    v1 = [0, v0, 0]
    v2 = 0.8*v1
    v3 = 1.2*v1

    pos = [r -r r]
    vel = hcat(v1, v2, v3)
    m = 0

    N = size(pos, 1)

    snap = lguys.Snapshot(pos, vel, zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
