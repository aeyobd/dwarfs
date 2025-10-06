using LilGuys
using PyFITS
import TOML
using Printf

function print_dm_props(galaxyname, orbitname, modelname)
    modeldir = joinpath(ENV["DWARFS_ROOT"], "analysis", galaxyname, orbitname, modelname)

    orb_props = TOML.parsefile(joinpath(modeldir, "orbital_properties.toml"))
    idx_f = orb_props["idx_f"]

    dm_scalars = read_fits(joinpath(modeldir, "profiles_scalars.fits"))
    @printf "r_circ_f_i\t%8.2f\n" dm_scalars.r_circ_max[idx_f]  / dm_scalars.r_circ_max[1]
    @printf "v_circ_f_i\t%8.2f\n" dm_scalars.v_circ_max[idx_f]  / dm_scalars.v_circ_max[1] 
    @printf "boundmass_f_i\t%8.6f\n" dm_scalars.bound_mass[idx_f] / dm_scalars.bound_mass[1] 
    #@printf "r_circ_i\t%8.2f\n" dm_scalars.r_circ_max[1] 
    #@printf "r_circ_f\t%8.2f\n" dm_scalars.r_circ_max[idx_f] 
    #@printf "v_circ_i\t%8.2f\n" dm_scalars.v_circ_max[1] * V2KMS
    #@printf "v_circ_f\t%8.2f\n" dm_scalars.v_circ_max[idx_f] * V2KMS
    #@printf "boundmass_i\t%8.6f\n" dm_scalars.bound_mass[1] 
    #@printf "boundmass_f\t%8.6f\n" dm_scalars.bound_mass[idx_f] 
end


println("scl smallperi")
print_dm_props("sculptor", "1e7_new_v31_r3.2", "orbit_smallperi")

println()
println("scl lmc-flyby")
print_dm_props("sculptor", "1e7_new_v25_r2.5", "smallperilmc")

println()
println("umi smallperi")
print_dm_props("ursa_minor", "1e7_new_v38_r4.0", "orbit_smallperi.5")
