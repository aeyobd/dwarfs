using LilGuys
using PyFITS
import TOML
using Printf

orb_props = TOML.parsefile("orbital_properties.toml")
idx_f = orb_props["idx_f"]

dm_scalars = read_fits("profiles_scalars.fits")
@printf "r_circ_i\t%8.2f\n" dm_scalars.r_circ_max[1] 
@printf "r_circ_f\t%8.2f\n" dm_scalars.r_circ_max[idx_f] 
@printf "v_circ_i\t%8.2f\n" dm_scalars.v_circ_max[1] * V2KMS
@printf "v_circ_f\t%8.2f\n" dm_scalars.v_circ_max[idx_f] * V2KMS
@printf "boundmass_i\t%8.6f\n" dm_scalars.bound_mass[1] 
@printf "boundmass_f\t%8.6f\n" dm_scalars.bound_mass[idx_f] 
