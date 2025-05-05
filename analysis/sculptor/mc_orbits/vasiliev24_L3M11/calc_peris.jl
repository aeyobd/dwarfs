using LilGuys
using PyFITS
using DataFrames, CSV

filename = "combined.hdf5"
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
