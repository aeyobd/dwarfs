using LilGuys
using PyFITS
import Agama
using Printf
import TOML


galaxyname = "ursa_minor"

df = read_fits("orbital_properties.fits")

pot = Agama.Potential(file = "agama_potential.ini")

function calc_r_J(halo, rho_peri)
    return LilGuys.find_zero(r -> LilGuys.mean_density(halo, r) - 3*rho_peri, halo.r_s)
end

halo = LilGuys.NFW(v_circ_max=31 / V2KMS, r_circ_max=3.2)


rho_peri = Agama.enclosed_mass(pot, df.pericentre) ./ (4π/3 * df.pericentre .^ 3)

r_J = calc_r_J.(halo, rho_peri)
r_J_arcmin = LilGuys.kpc2arcmin.(r_J, df.distance)


@printf "r_J kpc\t%0.2f ± %0.2f\n" LilGuys.mean(r_J) LilGuys.std(r_J)
@printf "r_J arcmin\t%0.2f ± %0.2f\n" LilGuys.mean(r_J_arcmin) LilGuys.std(r_J_arcmin)


obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))

σv = obs_props["sigma_v"] .+ randn(size(df, 1)) * obs_props["sigma_v_err"]

r_break = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri)
r_break_arcmin = LilGuys.kpc2arcmin.(r_break, df.distance)

@printf "r_break kpc\t%0.2f ± %0.2f\n" LilGuys.mean(r_break) LilGuys.std(r_break)
@printf "r_break arcmin\t%0.2f ± %0.2f\n" LilGuys.mean(r_break_arcmin) LilGuys.std(r_break_arcmin)

