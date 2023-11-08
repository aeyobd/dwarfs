import numpy as np
from lilguys import coords, Snapshot, units

fornax_obs = coords.observation(
    ra = 39.9970833,
    dec = -34.4491667,
    distance = 147,
    distance_err = 12,
    pm_ra = 0.374 ,
    pm_ra_err = 0.035 ,
    pm_dec = -0.401,
    pm_dec_err = 0.035,
    radial_velocity = 55.3,
    radial_velocity_err = 0.3
)


N = 10_000

mc_obs = coords.rand_coords(fornax_obs, N)
mc_phase = [coords.to_galcen(o) for o in mc_obs]
pos = np.array([[p.x, p.y, p.z] for p in mc_phase])
vel = np.array([[p.v_x, p.v_y, p.v_z] for p in mc_phase])
vel /= units.V_0
vel *= -1
pos /= units.R_0
snap = Snapshot(pos, vel, m=0)

snap.save("positions.hdf5")
