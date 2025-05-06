using LilGuys
using DataFrames

filename = "combined.hdf5"
out = Output(filename)

df = LilGuys.peris_apos(out)

LilGuys.write_fits("peris_apos.fits", df)
