using LilGuys
using DataFrames
using PyFITS

if isfile("peris_apos.fits")
    rm("peris_apos.fits")
end

filename = "combined.hdf5"
out = Output(filename)

df = LilGuys.peris_apos(out, verbose=true)

write_fits("peris_apos.fits", df)
