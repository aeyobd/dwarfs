import LilGuys as lguys


function main()
    r0 = 0.01
    v0 = 0.8

    pos = r0 .* [
        1  0   0  1 -1
        0  1   0  1 -0.5
        0  0  -1  1  0
        ]


    m = 0
    N = size(pos, 1)
    vel = zeros(3, N)

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
