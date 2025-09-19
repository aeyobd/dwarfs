using LilGuys
using CSV
using HDF5

module AgamaUtils
    include(joinpath(ENV["DWARFS_ROOT"], "utils/agama_utils.jl"))
end


let

    local times_code
    h5open("centres.hdf5", "r") do cen_file
        times_code = LilGuys.get_vector(cen_file, "times")
    end

    times_v21 = times_code * T2GYR / AgamaUtils.V_T2GYR
    println(times_v21)

    open("agama_times.txt", "w") do f
        for time in times_v21
            println(f, time)
        end
    end
end
