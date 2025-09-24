using LilGuys
using PyFITS
import Agama
using Printf
import TOML


function calc_r_J(halo, rho_peri)
    return 10 ^ LilGuys.find_zero(lr -> LilGuys.mean_density(halo, 10^lr) - 3*rho_peri, log10(halo.r_s))
end

function mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end

function print_radius(label, r, df)
    r_arcmin = LilGuys.kpc2arcmin.(r, df.distance)
    @printf "%s kpc   \t%0.2f ± %0.2f\n" label LilGuys.mean(r) LilGuys.std(r)
    @printf "%s arcmin\t%0.2f ± %0.2f\n" label LilGuys.mean(r_arcmin) LilGuys.std(r_arcmin)
end



function (@main)(ARGS)
    galaxyname = "ursa_minor"

    halo = LilGuys.NFW(v_circ_max=38 / V2KMS, r_circ_max=4.0)
    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))


    units = Agama.VASILIEV_UNITS
    df = read_fits("orbital_properties.fits")

    pot = Agama.Potential(file = "potential_mw_init.ini")
    pot_lmc = Agama.Potential(file = "potential_lmc_init.ini")



    rho_peri = mean_density(pot, df.pericentre, units)
    rho_peri_lmc = mean_density(pot_lmc, df.pericentre_lmc, units)

    r_J = calc_r_J.(halo, rho_peri)
    r_J_lmc = calc_r_J.(halo, rho_peri_lmc)
    print_radius("r_J", r_J, df)
    print_radius("r_J_lmc", r_J_lmc, df)


    σv = obs_props["sigma_v"] .+ randn(size(df, 1)) * obs_props["sigma_v_err"]

    r_break = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri)
    r_break_lmc = LilGuys.break_radius.(σv / V2KMS, df.time_last_peri_lmc)
    print_radius("r_break", r_break, df)
    print_radius("r_break_lmc", r_break_lmc, df)
end
