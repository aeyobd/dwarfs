import LilGuys as lguys
using Distributions 
using DataFrames
import TOML
import Random


mean_frame = lguys.GalactocentricFrame()

frame_distributions = Dict{Symbol, Distribution}(
    :d => Normal(mean_frame.d, 0.033),
    :ra => Normal(mean_frame.ra, 0.0001),
    :dec => Normal(mean_frame.dec, 0.00001),
    :z_sun => Normal(mean_frame.z_sun, 0.0003),
    :v_sun_x => Normal(mean_frame.v_sun[1], 3.0),
    :v_sun_y => Normal(mean_frame.v_sun[2], 1.4),
    :v_sun_z => Normal(mean_frame.v_sun[3], 0.08),
   )



function sample(obs_props::AbstractDict, frames_df::AbstractDict{<:Symbol, <:Distribution}, N::Integer=10_000)
    mc_obs = lguys.rand_coords(obs_props, N)
    frames = Vector{lguys.GalactocentricFrame}(undef, N)
    frames_df = DataFrame(Dict(
        :d => Vector{Float64}(undef, N),
        :ra => Vector{Float64}(undef, N),
        :dec => Vector{Float64}(undef, N),
        :z_sun => Vector{Float64}(undef, N),
        :v_sun_x => Vector{Float64}(undef, N),
        :v_sun_y => Vector{Float64}(undef, N),
        :v_sun_z => Vector{Float64}(undef, N)
       ))


    for i in 1:N
        for key in keys(frame_distributions)
            frames_df[i, key] = rand(frame_distributions[key])
        end

        frames[i] = lguys.GalactocentricFrame(
            d = frames_df[i, :d],
            ra = frames_df[i, :ra],
            dec = frames_df[i, :dec],
            z_sun = frames_df[i, :z_sun],
            v_sun = [
                frames_df[i, :v_sun_x],
                frames_df[i, :v_sun_y],
                frames_df[i, :v_sun_z]
            ]
        )
    end

    icrs_df = LilGuys.to_frame(mc_obs)

    # transform to phase space coordinates
    coords = [lguys.transform(lguys.Galactocentric, o, frame=frame) for (o, frame) in zip(mc_obs, frames)]
    return coords, frames_df, icrs_df
end


include("../../orbit_utils.jl")

function (@main)(ARGS)
    Random.seed!(496)

    units = get_units(".")
    pot = get_potential(".")
    obs_props = get_obs_props(".", "bootes3")
    coords, df_frames, df_icrs = sample(obs_props, frame_distributions, 1_000)

    orbits = LilGuys.agama_orbit(pot, coords; timerange=(0, -10/T2GYR), N=1000, agama_units=units)
    df_props = orbital_properties(pot, orbits, agama_units=units)

    df_galcen = LilGuys.to_frame(coords)

    df_all = hcat(df_props, df_galcen, df_frames)

    @info "writing properties"
    write_fits("orbital_properties.fits", df_all)
    @info "writing orbits"
    write_orbits(".", orbits)
end


