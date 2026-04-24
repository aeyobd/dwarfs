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


function get_samples(galaxyname, filename="mcmc_2exp")
    props = CSV.read( joinpath(ENV["DWARFS_ROOT"], "observations", 
                                galaxyname, "mcmc", "samples.$filename.csv"),
                    DataFrame)
end

function get_r_trans(prof, prof_outer)
	r_trans = LilGuys.find_zero(r -> LilGuys.surface_density(prof, r) - LilGuys.surface_density(prof_outer, r), 10)
    return r_trans
end
	


function get_r_trans(df)
	N = size(df, 1)
	r_trans = zeros(N)

	profs, prof_outer = get_profiles(df)
    return get_r_trans.(profs, prof_outer)
end

function get_profiles(df, dist=nothing)
    R_1 = df[!, "R_s"]
    f_outer = df[!, "f_outer"]
    R_2 =  df[!, "R_s_outer"]

	if !isnothing(dist)
		R_1 = LilGuys.arcmin2kpc(R_1, dist)
		R_2 = LilGuys.arcmin2kpc(R_2, dist)
	end
	
    prof = [LilGuys.Exp2D(R_s=R_1[i], M=(1-f_outer[i])) for i in eachindex(R_1)]
    prof_outer = [LilGuys.Exp2D(R_s=R_2[i], M=f_outer[i]) for i in eachindex(R_1)]

	return prof, prof_outer
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
    samples = get_samples(galaxyname)
    r_trans = get_r_trans(samples)
    print_quantity("R_trans", r_trans)
    print_param(df, "ellipticity_outer")
    print_param(df, "position_angle_outer")

    print_param(df_fe, "mu_fe_b")
    print_param(df_fe, "sigma_fe_b")
    println()
    println()
end


print_galaxy("sculptor")
print_galaxy("ursa_minor")
