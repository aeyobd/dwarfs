import LilGuys as lguys


function main()
    a = 3.944
    b = 0.311
    M = 6
    r_0 = a + b

    N = 1000
    v_circ(r) = sqrt(lguys.G * M * r^2 / (r^2 + r_0^2)^(3/2))

    r0 = 1
    δv0 = 0.01
    δr0 = 0.01

    rs =  r0 .+ r0/5 .* randn(N)
    vs = v_circ.(rs)

    pos = reshape(rs, 1, N) .* [1, 0, 0]
    vel = reshape(vs, 1, N) .* [0, 1, 0]

    pos .+= δr0 .* randn(3, N)
    vel .+= δv0 .* randn(3, N)

    m = 0

    snap = lguys.Snapshot(positions=pos, velocities=vel, masses=zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
