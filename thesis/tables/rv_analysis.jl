using LilGuys
import TOML
import StatsBase: quantile

module RVUtils
    include(joinpath(ENV["DWARFS_ROOT"], "observations/rv_utils.jl"))
end

function get_summary(galaxyname, modelname)
    df = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "velocities/processed/", "mcmc_properties_both.$modelname.toml"))
end

function rv_correction(galaxyname)
    obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxyname, "observed_properties.toml"))
    RVUtils.rv_gsr_shift(obs_props["ra"], obs_props["dec"])
end


function print_column(galaxyname, modelname)
    summary = get_summary(galaxyname, modelname)
    println(galaxyname, ":\t", modelname)

    df = summary["summary"]

    for param in ["μ", "σ", "dlσ_dlR"]
        i = findfirst(df["parameters"] .== [param])
        m = df["median"][i]
        if param == "μ"
            m += rv_correction(galaxyname)
        end

        print_quantity(param, m, df["error_lower"][i], df["error_upper"][i])
    end

    print_quantity("R grad", summary["R_grad_median"], summary["R_grad_el"], summary["R_grad_ep"])
    print_quantity("theta grad", summary["theta_grad_median"], summary["theta_grad_el"], summary["theta_grad_ep"])

    println("        bf sigma\t", round(summary["bf_rell"], digits=1))
    println("         bf grad\t", round(summary["bf_gradient"], digits=1))
    println()
end

function print_quantity(key, m, l, u)
    m = round(m, sigdigits=5)
    l = round(l, sigdigits=4)
    u = round(u, sigdigits=4)
    key = lpad(key, 16)
    println("$key\t$m  - $l  + $u")
end


print_column("sculptor", "rv_combined_x_wide_2c_psat_0.2")
print_column("sculptor", "rv_combined_x_wide_2c_psat_0.2_bin")
print_column("sculptor", "rv_tolstoy+23_x_wide_2c_psat_0.2")
print_column("sculptor", "rv_walker+09_x_wide_2c_psat_0.2")
print_column("sculptor", "rv_apogee_x_wide_2c_psat_0.2")

println()
print_column("ursa_minor", "rv_combined_x_2c_psat_0.2")
print_column("ursa_minor", "rv_combined_x_2c_psat_0.2_bin")
print_column("ursa_minor", "rv_pace+20_x_2c_psat_0.2")
print_column("ursa_minor", "rv_spencer+18_x_2c_psat_0.2")
print_column("ursa_minor", "rv_apogee_x_2c_psat_0.2")
