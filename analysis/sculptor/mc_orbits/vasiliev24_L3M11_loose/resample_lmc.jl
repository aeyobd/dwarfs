using CSV, DataFrames
using LilGuys


lmc_file = ENV["DWARFS_ROOT"] * "/simulations/sculptor/mc_orbits/vasiliev_lmc/trajlmc.txt"

lmc_traj = CSV.read(lmc_file, DataFrame, delim=" ", header = [:time, :x, :y, :z, :v_x, :v_y, :v_z])



V_T2GYR = 0.97779

# times are run in reverse for gadget's sake
lmc_traj.time = -lmc_traj.time * V_T2GYR / T2GYR


out = Output(".")
println(length(out.times))

times = out.times


lmc_x = LilGuys.lerp(lmc_traj.time, lmc_traj.x).(times)
lmc_y = LilGuys.lerp(lmc_traj.time, lmc_traj.y).(times)
lmc_z = LilGuys.lerp(lmc_traj.time, lmc_traj.z).(times)

lmc_vx = LilGuys.lerp(lmc_traj.time, lmc_traj.v_x).(times)
lmc_vy = LilGuys.lerp(lmc_traj.time, lmc_traj.v_y).(times)
lmc_vz = LilGuys.lerp(lmc_traj.time, lmc_traj.v_z).(times)

lmc = DataFrame(
    :time => times,
    :x => lmc_x,
    :y => lmc_y,
    :z => lmc_z,
    :v_x => lmc_vx / V2KMS,
    :v_y => lmc_vy / V2KMS,
    :v_z => lmc_vz / V2KMS,
   )

CSV.write("lmc_traj.csv", lmc)
