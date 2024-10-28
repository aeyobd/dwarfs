using LilGuys
using DataFrames

filename = "combined.hdf5"
out = Output(filename)

idx, peris, apos = LilGuys.peris_apos(out)

snap1 = out[1]

df = DataFrame(index=idx, pericenter=peris, apocenter = apos)

idx_sort = sortperm(df.index)
df = df[idx_sort, :]

LilGuys.write_fits("peris_apos.fits", df)
