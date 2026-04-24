using LilGuys
import TOML
using PyFITS
import StatsBase: quantile
import DataFrames: disallowmissing
using Printf
import Agama
import HDF5

include("table_utils.jl")


function get_orbit_props(galaxyname, modelname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "simulations", galaxyname, modelname, "orbit.toml"))
    orbit_props
end


function get_orbit(galaxyname, modelname)
    orbit = LilGuys.Orbit(joinpath(ENV["DWARFS_ROOT"], "simulations", galaxyname, modelname, "orbit.csv"))
    orbit
end


function print_orbit(galaxyname, modelname, starsname)
    orbit_props = get_orbit_props(galaxyname, modelname)
    orbit = get_orbit(galaxyname, modelname)

    df_result = TOML.parsefile(
        joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname, "orbital_properties.toml")
       )

    @printf "t_i\t%0.2f\n" -get(df_result, "t_f_gyr", 0)
    @printf "x_i\t%0.2f\n" orbit.positions[1, 1]
    @printf "y_i\t%0.2f\n" orbit.positions[2, 1]
    @printf "z_i\t%0.2f\n" orbit.positions[3, 1]


    @printf "vx_i\t%0.2f\n" orbit.velocities[1, 1]*V2KMS
    @printf "vy_i\t%0.2f\n" orbit.velocities[2, 1]*V2KMS
    @printf "vz_i\t%0.2f\n" orbit.velocities[3, 1]*V2KMS

    @printf "pericentre\t%0.2f\n" LilGuys.pericenter(orbit)

    halo = LilGuys.load_profile(joinpath(ENV["DWARFS_SIMS"], galaxyname, modelname, "../halo.toml"))

    @printf "vmax\t%0.2f\n" LilGuys.v_circ_max(halo)*V2KMS
    @printf "rmax\t%0.2f\n" LilGuys.r_circ_max(halo)

    HDF5.h5open(joinpath(ENV["DWARFS_SIMS"], galaxyname, modelname, "out/snapshot_001.hdf5"), "r") do f
        header = HDF5.attrs(f["Header"])
        @printf "part mass\t%0.2e" header["MassTable"][2] * M2MSUN
    end

    println()
end




println("sculptor smallperi ")
print_orbit(modelnames["scl_smallperi"]...,)


println("scl lmc ")
print_orbit(modelnames["scl_lmc"]...,)

println("umi.5 ")
print_orbit(modelnames["umi_smallperi"]...,)

println("umi.1")
print_orbit("ursa_minor", "1e5_v38_r4.0/orbit_smallperi.1", "")
