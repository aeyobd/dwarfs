import numpy as np
from lilguys import coords, Snapshot, units

observations = coords.observation(
    ra = 15.03917,
    dec = -33.70917,
    distance = 86,
    distance_err = 6,
    pm_ra = 0.099,
    pm_ra_err = 0.002 ,
    pm_dec = -0.160,
    pm_dec_err = 0.002,
    radial_velocity = 111.4,
    radial_velocity_err = 0.37
)


N = 10_000

mc_obs = coords.rand_coords(observations, N)
mc_phase = [coords.to_galcen(o) for o in mc_obs]
pos = np.array([[p.x, p.y, p.z] for p in mc_phase])
vel = np.array([[p.v_x, p.v_y, p.v_z] for p in mc_phase])
vel *= -1
vel /= units.V_0

pos /= units.R_0
snap = Snapshot(pos, vel, m=0)

snap.save("positions.hdf5")
