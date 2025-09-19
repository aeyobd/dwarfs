using LilGuys

out = Output(ENV["DWARFS_ROOT"] * "/analysis/sculptor/1e6_V31_r3.2/vasiliev24_L3M11_extremeperi/")
snap = out[end]

snap.positions .-= snap.x_cen
snap.velocities .-= snap.v_cen

LilGuys.save("initial.hdf5", snap)
