using LilGuys
import StatsBase: sample
using HDF5

filename = "combined.hdf5"
out = Output(filename)

exclude_idx = 1

idx = sample(out[1].index[1+exclude_idx:end], 1000)


@info "Extracting trajectory"
positions = LilGuys.extract_vector(out, :positions, idx)
velocities = LilGuys.extract_vector(out, :velocities, idx)
times = out.times

filename = "trajectory.hdf5"

h5open(filename, "w") do file
    file["positions"] = positions
    file["velocities"] = velocities
    file["times"] = times
end