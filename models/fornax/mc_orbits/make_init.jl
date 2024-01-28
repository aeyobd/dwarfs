import LilGuys as lguys


function main()
    N = 1000

    r0 = 100
    v0 = 1
    σ = 0.05
    rs = r0 .+ σ*r0 .* randn(N)
    vs = v0 .+ σ*v0 .* randn(N)

    pos = reshape(rs, 1, N) .* lguys.rand_unit()
    vel = reshape(vs, 1, N) .* lguys.rand_unit()

    snap = lguys.Snapshot(pos, vel, zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
