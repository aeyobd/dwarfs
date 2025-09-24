using LilGuys
using PyFITS
import TOML
using Printf

orb_props = TOML.parsefile("../orbital_properties.toml")
idx_f = orb_props["idx_f"]

filename = joinpath(ARGS[1], "stellar_profiles_3d_scalars.fits")
df = read_fits(filename)

idx = argmin(abs.(df.snapshot .- idx_f))
@info "reading in $(df.snapshot[idx]), closest to $idx_f"

@printf "boundmass_i\t%8.6f\n" df.bound_mass[1] 
@printf "boundmass_f\t%8.6f\n" df.bound_mass[idx] 
@printf "sigma_v_i\t%8.6f\n" df.sigma_v[1] 
@printf "sigma_v_f\t%8.6f\n" df.sigma_v[idx] 
