#!/bin/env julia


using LilGuys

function main()
    out = Output("../out")

    pos = out.x_cen
    vel = out.v_cen
    
    snap = Snapshot(pos, vel, zeros(size(pos, 2)))
    df = LilGuys.to_gaia(snap, p_min=0, SkyFrame=LilGuys.GSR)

    LilGuys.write_fits("sky_orbit.fits", df)

end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
