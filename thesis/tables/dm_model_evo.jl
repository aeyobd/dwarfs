using PyFITS
import TOML
using Printf
include("utils.jl")

function print_dm_props(galaxyname, modelname, starsname)
    modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, modelname)

    orb_props = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))
    idx_f = orb_props["idx_f"]

    dm_scalars = read_fits(joinpath(modeldir, "profiles_scalars.fits"))
    @printf "r_circ_f_i\t%8.3f\n" dm_scalars.r_circ_max[idx_f]  / dm_scalars.r_circ_max[1]
    @printf "v_circ_f_i\t%8.3f\n" dm_scalars.v_circ_max[idx_f]  / dm_scalars.v_circ_max[1] 
    @printf "boundmass_f_i\t%8.6f\n" dm_scalars.bound_mass[idx_f] / dm_scalars.bound_mass[1] 
end


println("scl smallperi")
print_dm_props(modelnames["scl_smallperi"]...)

println()
println("scl lmc-flyby")
print_dm_props(modelnames["scl_lmc"]...)

println()
println("umi smallperi")
print_dm_props(modelnames["umi_smallperi"]...)
