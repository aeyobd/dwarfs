using LilGuys
import Random
Random.seed!(496)


import TOML


function get_obs_props(galaxyname)
    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
end

function sample(galaxyname)
    obs_props = get_obs_props(galaxyname)
    N = 16000

    Mv = obs_props["Mv"] .+ randn(N) * obs_props["Mv_err"]
    Lv = LilGuys.mag_to_L.(Mv)
    log_M_L_err = obs_props["M_L_s_err"] / obs_props["M_L_s"] / log(10)
    M_L = obs_props["M_L_s"]  * (10 .^ (log_M_L_err*randn(N)))
    Mstar = Lv .* M_L



    vmax = vel_from_M_s.(Mstar / M2MSUN) .* 10 .^ (0 .+ 0.035 * randn(N))
    rmax = LilGuys.Ludlow.solve_rmax.(vmax, 0.1*randn(N),)

    halos = [LilGuys.NFW(r_circ_max=r, v_circ_max=v) for (r,v) in zip(rmax, vmax)]


    print_quantity("Lv", Lv)
    print_quantity("M_L", M_L)
    print_quantity("Mstar", Mstar)
    print_quantity("vmax", vmax * V2KMS)
    print_quantity("rmax", rmax)
    print_quantity("M200", LilGuys.M200.(halos))
    print_quantity("c", [h.c for h in halos])
end

function print_quantity(label, quantity)
    l, m, h = LilGuys.quantile(quantity, [0.16, 0.5, 0.84])
    m = LilGuys.Measurement(m, m-l, h-m)
    println(label, "\t", m)
end

vel_from_M_s(ms) = LilGuys.find_zero(x->LilGuys.M_s_from_vel_fattahi(x) - ms, 0.1)

println("sculptor")
sample("sculptor")
println()
println("ursa_minor")
sample("ursa_minor")
