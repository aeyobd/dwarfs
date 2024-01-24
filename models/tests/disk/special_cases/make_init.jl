import LilGuys as lguys


function main()
    r0 = 8
    v0 = 0.8

    pos = r0 .* [
        1  1.2    0  
        0  0      0 
        0  0.05   1
        ]


    vel = v0 .* [
        0    0    0
        1    1.2  0
        0.05 0.05  0
        ]

    m = 0
    N = size(pos, 1)

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
