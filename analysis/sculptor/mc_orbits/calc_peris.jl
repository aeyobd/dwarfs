using LilGuys
using DataFrames
using PyFITS
rm("peris_apos.fits")

filename = "combined.hdf5"
out = Output(filename)

df = LilGuys.peris_apos(out, verbose=true)

write_fits("peris_apos.fits", df)
