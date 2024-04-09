using LilGuys
using CSV
using DataFrames

filename = "out"
out = Output(filename)

idx, peris, apos = LilGuys.peris_apos(out)

snap1 = out[1]

df = DataFrame(index=idx, pericenter=peris, apocenter = apos)
CSV.write("peris_apos.csv", df)

println("finished")

