import LilGuys as lguys


function main()
    p0 = 1
    v0 = 0.3
    r1 = [0.97000436, -0.24308753, 0]
    r2 = -r1
    r3 = zeros(3)

    v3 = [-0.93240737, -0.86473146, 0]
    v1 = -v3/2
    v2 = -v3/2

    pos = hcat(r1, r2, r3)
    vel = hcat(v1, v2, v3)
    m = 1

    s = lguys.Snapshot(pos, vel, m)
    lguys.save("threebody.hdf5", s)
end


main()
