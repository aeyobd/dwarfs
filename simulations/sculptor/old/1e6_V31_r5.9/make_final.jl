using LilGuys

out = Output(".")

snap_f = out[end]

LilGuys.save("final.hdf5", snap_f, centre=true)

snap_mid = out[20]
LilGuys.save("mid.hdf5", snap_mid, centre=true)
