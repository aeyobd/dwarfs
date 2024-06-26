import LilGuys as lguys

function main()
    p0 = 1
    v0 = 0.3
    x0 = [p0 -p0
          0   0 
          0   0] 

    v0 = [0 0
          v0 -v0
          0 0]

    m = [1, 1]

    s = lguys.Snapshot(x0, v0, m)
    lguys.save("twobody.hdf5", s)
end


main()
