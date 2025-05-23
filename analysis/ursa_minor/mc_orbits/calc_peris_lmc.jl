using LilGuys
using DataFrames, CSV
using PyFITS



if length(ARGS) > 0
    t_max = ARGS[1]
else
    t_max = 0
end

filename = "combined.hdf5"
if isfile("peris_apos.fits")
    rm("peris_apos.fits")
end

out = Output(filename)

traj_lmc = CSV.read("lmc_traj.csv", DataFrame)
pos_lmc = [traj_lmc.x traj_lmc.y traj_lmc.z]'


df = LilGuys.peris_apos(out)

df_lmc = LilGuys.peris_apos(out, x0 = pos_lmc)

df[!, :peri_lmc] = df_lmc.pericentre
df[!, :apo_lmc] = df_lmc.apocentre
df[!, :t_last_peri_lmc] = df_lmc.t_last_peri
df[!, :t_last_apo_lmc] = df_lmc.t_last_apo

write_fits("peris_apos.fits", df)
