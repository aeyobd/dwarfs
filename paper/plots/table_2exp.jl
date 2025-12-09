using LilGuys
import TOML
using Printf
import Agama
using DataFrames, CSV

include("table_utils.jl")

α = LilGuys.R_h(LilGuys.Exp2D())

function get_summary(galaxyname, filename="mcmc_2exp")
    props = CSV.read( joinpath(ENV["DWARFS_ROOT"], "observations", 
                                galaxyname, "mcmc", "summary.$filename.csv"),
                    DataFrame)
end

function get_metallicity_fit(galaxyname)
    return CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/multipop/processed", "$galaxyname.mcmc_2pop_vel_fe.summary.csv"),
                    DataFrame)
end

function get_param(df, key)
    idx = only(findall(df.parameters .== [key]))

    return Measurement(df.median[idx], df.lower_error[idx], df.upper_error[idx])
end

function print_param(df, key)
    print_quantity(rpad(key, 20), get_param(df, key))
end

function print_R_h(df, key)
    print_quantity(replace(rpad(key, 20), "R_s"=>"R_h"), get_param(df, key) * α)
end

function print_galaxy(galaxyname)
    println(galaxyname)
    println("-"^30)
    df = get_summary(galaxyname)
    df_fe = get_metallicity_fit(galaxyname)
    print_param(df, "R_s")
    print_R_h(df, "R_s")
    print_param(df, "ellipticity")
    print_param(df, "position_angle")
    print_param(df_fe, "mu_fe_a")
    print_param(df_fe, "sigma_fe_a")
    println()
    print_param(df, "f_outer")
    print_param(df, "R_s_outer")
    print_R_h(df, "R_s_outer")
    print_param(df, "ellipticity_outer")
    print_param(df, "position_angle_outer")

    print_param(df_fe, "mu_fe_b")
    print_param(df_fe, "sigma_fe_b")
    println()
    println()
end


print_galaxy("sculptor")
print_galaxy("ursa_minor")
