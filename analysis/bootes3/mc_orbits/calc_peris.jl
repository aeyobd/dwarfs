using LilGuys
using DataFrames
using PyFITS

filename = "combined.hdf5"
out = Output(filename)

df = LilGuys.peris_apos(out)

write_fits("peris_apos.fits", df)
