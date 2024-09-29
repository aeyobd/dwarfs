import LilGuys as lguys


function main()
    N = 1000

    r0 = 100
    v0 = 1
    rs = r0 .+ r0/5 .* randn(N)
    vs = v0 .+ v0/5 .* randn(N)

    pos = reshape(rs, 1, N) .* lguys.rand_unit(N)
    vel = reshape(vs, 1, N) .* lguys.rand_unit(N)

    snap = lguys.Snapshot(pos, vel, zeros(N))
    lguys.save("initial.hdf5", snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
