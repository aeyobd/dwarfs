using LilGuys

using CSV, DataFrames
import TOML
using HDF5

using PythonCall
agama = pyimport("agama")
np = pyimport("numpy")


include(ENV["DWARFS_ROOT"] * "/utils/agama_utils.jl")
include(ENV["DWARFS_ROOT"] * "/utils/agama_plots.jl")
include(ENV["DWARFS_ROOT"] * "/utils/lmc_utils.jl")


function scl_lmc_sample_and_orbit(Φ, scl, scl_err, lmc, lmc_err; 
        Mlmc=20, N=1000, units=:vasiliev, 
        Mlmc_err=0.2, a_fric=1.1, a_fric_err=0.2, dt_max=5, 
        t_max = -1 / T2GYR,
        kwargs...
    )
    sigmafnc = calc_σv_func(Φ)
    
    icrs_samples = LilGuys.rand_coords(lmc, lmc_err, N)
    gc_samples = LilGuys.transform.(Galactocentric, icrs_samples)
    icrs_scl_samples = LilGuys.rand_coords(scl, scl_err, N);
    gc_scl_samples = LilGuys.transform.(Galactocentric, icrs_scl_samples);
    Mlmcs = Mlmc * 10 .^ (Mlmc_err * randn(N))
    r_s = @. 8.5 * (Mlmc / 10) ^ 0.6

    icrs_df = LilGuys.to_frame(icrs_samples)
    icrs_scl_df = LilGuys.to_frame(icrs_scl_samples)
    a_fric_samples = a_fric .* 10 .^ (a_fric_err .* randn(N))


    orbits = Orbit[]
    orbits_scl = Orbit[]
    Φs = []
    if units == :vasiliev
        time_scale = T2GYR / V_T2GYR
    else
        time_scale = 1
    end
    
    for i in 1:N
        print("orbit $i\r")
        gc = gc_samples[i]
        m = Mlmcs[i]
        
        Φ_new, orbit = make_lmc_mw_pot(Φ, gc, 
            time=t_max, 
            units = units, 
            Mlmc=t->m,  
            σv=x->py2f(sigmafnc(x)), 
            r_s=r_s, 
            reflex_motion=true, 
            dynamical_friction=a_fric_samples[i], 
            dt_max=dt_max
        ) 
        
        push!(orbits, orbit)
        push!(Φs, Φ_new)

        
        gc = gc_scl_samples[i]
        orbit = calc_orbit(gc, Φs[i], time=orbits[i].time * time_scale , units=units);
        push!(orbits_scl, orbit)
    end

    orbits_scl_lmc = [orbits_scl[i] - orbits[i] for i in eachindex(orbits)]


    # properties
    peri_lmc = Vector{Float64}(undef, N)
    peris = Vector{Float64}(undef, N)
    t_peri = Vector{Float64}(undef, N)
    t_peri_lmc = Vector{Float64}(undef, N)
    
    for i in 1:N
        r = calc_r(orbits_scl[i].position)
        idx = argmin(r)
        t_peri[i] = orbits_scl[i].time[idx]
        peris[i] = r[idx]
    
        r = calc_r(orbits_scl_lmc[i].position)
        idx = argmin(r)
        t_peri_lmc[i] = orbits_scl_lmc[i].time[idx]
        peri_lmc[i] = r[idx]
    end

    # df all
    df_all = copy(icrs_scl_df)
    for symbol in names(df_all)
        df_all[:, "$(symbol)_lmc"] = icrs_df[:, symbol]
    end
    
    df_all[:, "t_peri"] = t_peri * T2GYR
    df_all[:, "t_peri_lmc"] = t_peri_lmc * T2GYR
    df_all[:, "peri"] = peris
    df_all[:, "peri_lmc"] = peri_lmc
    df_all[:, "Mlmc"] = Mlmcs
    df_all[:, "a_fric"] = a_fric_samples
        
    return orbits, orbits_scl, orbits_scl_lmc, Φs, df_all
end


function orbits_to_hdf5(filename, orbits)
    N = length(orbits)

    h5open(filename, "w") do file
        for i in 1:N
            orbit = orbits[i]
            group = string("orbit_", i)
            g = HDF5.create_group(file, group)
            g["time"] = orbit.time
            g["position"] = orbit.position
            g["velocity"] = orbit.velocity
        end
    end
end


function main()
    if length(ARGS) < 1
        println("Usage: julia calc_orbits.jl <dirname>")
        return
    end

    dir = ARGS[1]
    cd(dir)

    Φ = agama.Potential("agama_potential.ini")
    icrs_scl = LilGuys.coord_from_file("scl.toml")
    icrs_scl_err = LilGuys.coord_err_from_file("scl.toml")
    icrs_lmc = LilGuys.coord_from_file("lmc.toml")
    icrs_lmc_err = LilGuys.coord_err_from_file("lmc.toml")

    params = TOML.parsefile("params.toml")
    N = get(params, "N", 100)
    units = get(params, "units", :code) |> Symbol
    dt_max = get(params, "dt_max", 100)
    err_scale = get(params, "err_scale", 1.0)
    Mlmc = get(params, "Mlmc", 20)
    Mlmc_err = get(params, "Mlmc_err", 0.2)
    a_fric = get(params, "a_fric", 1.1)
    a_fric_err = get(params, "a_fric_err", 0.2)
    reflex_motion = get(params, "reflex_motion", true)


    icrs_scl_err = ICRS(
        ra=icrs_scl_err.ra * err_scale, 
        dec=icrs_scl_err.dec * err_scale, 
        distance=icrs_scl_err.distance * err_scale,
        pmra=icrs_scl_err.pmra * err_scale,
        pmdec=icrs_scl_err.pmdec * err_scale,
        radial_velocity=icrs_scl_err.radial_velocity * err_scale
    )

    icrs_lmc_err = ICRS(
        ra=icrs_lmc_err.ra * err_scale, 
        dec=icrs_lmc_err.dec * err_scale, 
        distance=icrs_lmc_err.distance * err_scale,
        pmra=icrs_lmc_err.pmra * err_scale,
        pmdec=icrs_lmc_err.pmdec * err_scale,
        radial_velocity=icrs_lmc_err.radial_velocity * err_scale
    )

    orbits, orbits_scl, orbits_scl_lmc, Φs, df_all = scl_lmc_sample_and_orbit(
        Φ, 
        icrs_scl, icrs_scl_err, 
        icrs_lmc, icrs_lmc_err, 
        N=N, units=units, dt_max=dt_max,
        Mlmc=Mlmc, Mlmc_err=Mlmc_err, a_fric=a_fric, a_fric_err=a_fric_err,
        reflex_motion=reflex_motion,
    )

    CSV.write("perisapos.csv", df_all)
    orbits_to_hdf5("orbits_lmc.hdf5", orbits)
    orbits_to_hdf5("orbits_scl.hdf5", orbits_scl)
    orbits_to_hdf5("orbits_scl_lmc.hdf5", orbits_scl_lmc)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
