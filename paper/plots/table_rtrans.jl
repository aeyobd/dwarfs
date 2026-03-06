using LilGuys
import TOML
using Printf
import Agama
using DataFrames, CSV

include("table_utils.jl")

α = LilGuys.R_h(LilGuys.Exp2D())

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
    df = get_samples(galaxyname)
    r_trans = get_r_trans(df)
    print_quantity("r_trans", r_trans)
    println()
    println()
end


print_galaxy("sculptor")
print_galaxy("ursa_minor")
