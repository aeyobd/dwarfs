using PyFITS
using LilGuys

include("utils.jl")

function load_model(galaxyname, modelname)
    return read_fits(joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, modelname, "orbital_properties.fits"))
end


function print_period(galaxyname, modelname)
    df = load_model(galaxyname, modelname)
    println(galaxyname, " ", modelname)
    print_quantity("period", collect(skipmissing(df.period)) * T2GYR)
    println("number missing:\t", LilGuys.mean(ismissing.(df.period)))
    println()
end

print_period("sculptor", "EP2020")
print_period("sculptor", "vasiliev24_L3M11_9Gyr")
print_period("sculptor", "L3M10")
print_period("sculptor", "L2M11")
