using LilGuys
import TOML
using PyFITS
import StatsBase: quantile
import DataFrames: disallowmissing
using Printf
import Agama
import HDF5

include("table_utils.jl")


function get_stars_prof(galaxyname, modelname, starsname)
    stars_prof = LilGuys.load_profile(joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "../stars", starsname, "profile.toml"))
end


function print_orbit(galaxyname, modelname, starsname)
    println(galaxyname, modelname, starsname)
    stars_prof = get_stars_prof(galaxyname, modelname, starsname)
    @printf "R_h\t%4.3f\n" LilGuys.R_h(stars_prof)

    println()
end




println("sculptor smallperi ")
print_orbit(modelnames["scl_smallperi"]...,)
print_orbit(modelnames["scl_smallperi_plummer"]...,)


println("scl lmc ")
print_orbit(modelnames["scl_lmc"]...,)
print_orbit(modelnames["scl_lmc_plummer"]...,)

println("umi.5 ")
print_orbit(modelnames["umi_smallperi"]...,)
print_orbit(modelnames["umi_smallperi_plummer"]...,)

