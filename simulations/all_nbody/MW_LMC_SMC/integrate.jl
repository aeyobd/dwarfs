import Agama
import LilGuys
import Random
using LilGuys

using CSV, DataFrames

Random.seed!(127)
galaxynames = ["mw", "lmc", "smc", "sculptor", "ursa_minor"]


module NBody 
	include("../nbody_utils.jl")
    include("../initial_conditions.jl")
end


function get_obs_props()
    df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv"), DataFrame)

    df[!, :distance_modulus] = LilGuys.kpc2dm.(df.distance)
    df[!, :ra_err] .= 0
    df[!, :dec_err] .= 0
    return df
end


const OBS_PROPS = get_obs_props()

function get_obs_props(galaxyname)
    df = OBS_PROPS
    idx = findfirst(df.galaxyname .== [galaxyname])
    row = df[idx, :]
    return Dict((k=>row[k] for k in names(row)))
end

POT_STARS = NBody.get_potential("EP2020_stars")


function sample_posvel(obs_props)
    icrs = LilGuys.rand_coord(obs_props)
    gc_i = LilGuys.transform(Galactocentric, icrs)

    pos = NBody.Point(LilGuys.position(gc_i))
    vel = NBody.Point(LilGuys.velocity(gc_i) / V2KMS)
    return pos, vel
end

function sample_mw_halo()
    # correamagnus+vasiliev2022
    # todo: include halo slope uncertanties
    M200_mw = 10 ^ (2.04  + 0.10 * randn())
    r_s_mw = 10 ^ (1.11 + 0.5 * randn())

    halo_mw_dm = NFW(M200=M200_mw, r_s=r_s_mw)
    halo_mw_dm
end

function sample_lmc_halo()
    #koposov
    M_32_8_lmc = 5.75 + 0.9 * randn()
    r_s_lmc = 2.71 * 10^(0.29 * randn())
    halo_lmc_unscaled = NFW(M_s=1.0, r_s=r_s_lmc, c=1)

    M_s_lmc = M_32_8_lmc / mass(halo_lmc_unscaled, 32.8)
    halo_lmc = NFW(M_s=M_s_lmc, r_s=r_s_lmc)

    return halo_lmc
end

function sample_smc_halo()
    vmax_smc = (56 + 5* randn()) / V2KMS
    σ_ludlow = 0.1
    rmax_smc = LilGuys.Ludlow.solve_rmax(vmax_smc, σ_ludlow * randn())
    halo_smc = NFW(r_circ_max=rmax_smc, v_circ_max=vmax_smc)

    return halo_smc
end

function sample_halo(obs_props)
    σ_fattahi = 0.037
    σ_ludlow = 0.1
    σ_M_L = 0.2

    M_L_star = 2.0 * 10^(σ_M_L * rand())
    L = LilGuys.mag_to_L(obs_props["Mv"])
    Mstar = L / M2MSUN * M_L_star

    vcircmax = LilGuys.vel_from_M_s_fattahi(Mstar * 10^(σ_fattahi*randn()))


    rcircmax = LilGuys.Ludlow.solve_rmax(vcircmax, σ_ludlow * randn())
    halo = LilGuys.NFW(v_circ_max=vcircmax, r_circ_max=rcircmax)
end



function sample_ic(galaxyname::String)
    if galaxyname == "mw"
        pos, vel =  NBody.Point(0.,0,0), NBody.Point(0.,0,0)
        halo = sample_mw_halo()
    else
        obs_props = get_obs_props(galaxyname)
        pos, vel = sample_posvel(obs_props)
        if galaxyname == "lmc"
            halo = sample_lmc_halo()
        elseif galaxyname == "smc"
            halo = sample_smc_halo()
        else
            halo = sample_halo(obs_props)
        end
    end

    return pos, vel, halo
end


function initial_conditions()
    ics = sample_ic.(galaxynames)
    ic_pos = [ic[1] for ic in ics]
    ic_vel = [ic[2] for ic in ics]
    ic_halo = [ic[3] for ic in ics]
    return ic_pos, ic_vel, ic_halo
end


function clear_dir(dir)
    if !isdir(dir)
        mkdir(dir)
    end
    for file in readdir(dir)
        rm(joinpath(dir, file))
    end
end


function (@main)(ARGS)
    units = Agama.AgamaUnits()
    obs_props = NBody.get_obs_props()
    Norbits = 100
    timestep = 0.02

    clear_dir("out")
    clear_dir("out_nosmc")
    clear_dir("out_nofric")

    for i in 1:Norbits
        @info "integrating orbit $i"

        ic_pos, ic_vel, ic_halo = initial_conditions()
        NBody.write_initial_conditions("out/initial_$i.fits", [ic_pos, ic_vel, ic_halo])

        halo_mw_dm = ic_halo[1]
        #halo_mw = Agama.Potential(type="NFW", mass=halo_mw_dm.M_s, scaleRadius=halo_mw_dm.r_s)  + POT_STARS

        #ic_halo = [halo_mw; ic_halo[2:end]]

        σv = x -> LilGuys.v_circ(halo_mw_dm, radii(x)) / √3
        ρ = x -> LilGuys.density(halo_mw_dm, radii(x))


        f_fric = NBody.make_dyn_fric_models(halo_mw_dm, units, ic_halo[2:end], f_σ=σv, f_ρ=ρ)
        f_fric = [NBody.Massless(); f_fric]

        function force(accelerations, positions, velocities, halos, times)
            NBody.nbody_acceleration!(accelerations, positions, halos)
            NBody.force_dyn_friction!(accelerations, positions .- [positions[1]], velocities .- [velocities[1]], f_fric)
        end

        orbits = NBody.integrate_particles(
            ic_pos, ic_vel, ic_halo, force,
            timestep=timestep,
        )
        NBody.write_orbits("out/orbits_$i.hdf5", orbits)

        idxs = galaxynames .!= ["smc"]
        orbits_nosmc = NBody.integrate_particles(
                                                 ic_pos[idxs], ic_vel[idxs], ic_halo[idxs], force,
            timestep=timestep,
        )
        NBody.write_orbits("out_nosmc/orbits_$i.hdf5", orbits_nosmc)


        function force_no(accelerations, positions, velocities, halos, times)
            NBody.nbody_acceleration!(accelerations, positions, halos)
        end
        orbits_no = NBody.integrate_particles(
            ic_pos, ic_vel, ic_halo, force_no,
            timestep=timestep,
        )
        NBody.write_orbits("out_nofric/orbits_$i.hdf5", orbits_no)
    end
end
