using LilGuys
import TOML
using PyFITS
import StatsBase: quantile
import DataFrames: disallowmissing
using Printf
import Agama

galaxyname = "bootes3"
include("table_utils.jl")



function print_radius(label, r, distance)
    r_arcmin = LilGuys.kpc2arcmin.(r, distance)
    print_quantity(label * " / kpc", r)
    print_quantity(label * " / arcmin", r_arcmin)
end

function mean_density(pot, r, units)
    M =  Agama.enclosed_mass(pot, r, units) 
    V = (4π/3 * r .^ 3)

    return M ./ V
end

function get_orbit(modelname)
    orbit_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "orbits", galaxyname, "EP2020_special_cases", "$modelname.toml"))
    orbit_props
end

column_labels = ["pericentre", "apocentre", "distance", "pmra", "pmdec", "radial_velocity", "period"] 


function print_orbit(modelname)
    orbit_props = get_orbit(modelname)

    @printf "%16s\t" modelname
    for key in column_labels
        val = orbit_props[key]
        if val isa Real
            val = round(orbit_props[key], sigdigits=5)
        end
        if key == "period"
            val = abs(val) * T2GYR
        end

        @printf "%8.2f\t" val
    end


    println()
end




@printf "%16s\t" "modelname"

for key in column_labels
    @printf"%8s\t" key
end
println()

for modelname in ["orbit_peri_1.5", "orbit_peri_4", "orbit_mean", 
                  "orbit_peri_12", "orbit_peri_18", "orbit_peri_26"]

    print_orbit(modelname)
end
