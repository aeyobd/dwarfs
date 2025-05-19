import LilGuys as lguys
using Distributions 
using DataFrames
import TOML


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


"""
    sample(obs_props, N=10_000)

Given a dictionary of observed properties
(ra, dec, pmra, pmdec, distance_modulus, radial_velocity) and uncertanties
samples N random coordinate, with the velocities backwards in time
Returns a snapshot of massless particles.
"""
function sample(obs_props::AbstractDict, N::Integer = 10_000)
    mc_obs = lguys.rand_coords(obs_props, N)

    mc_phase = lguys.transform.(lguys.Galactocentric, mc_obs)
    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N))

    return snap
end



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

    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o, frame=frame) for (o, frame) in zip(mc_obs, frames)]
    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N))

    return snap, frames_df
end



"""
    frame_from_params(params)

Retrieve the galactocentric frame specification if 
there is a frame dictionary in the params
"""
function frame_from_params(params)
    if "frame" in keys(params)

        df = params["frame"]
        frame = lguys.Galactocentric(
            galcen_distance = df["galcen_distance"],
            galcen_v_sun = df["galcen_v_sun"],
            z_sun = df["z_sun"],
            roll = df["roll"],
            galcen_coord = lguys.ICRS(
                ra = df["galcen_coord"]["ra"],
                dec = df["galcen_coord"]["dec"],
                distance = df["galcen_coord"]["distance"],
                pmra = df["galcen_coord"]["pmra"],
                pmdec = df["galcen_coord"]["pmdec"],
                radial_velocity = df["galcen_coord"]["radial_velocity"],
            )
        )
    else
        frame = lguys.default_gc_frame
    end

    return frame
end


function sample_special_cases(filename="initial_conditions.toml")
    df = TOML.parsefile(filename)

    obs_special = lguys.ICRS[]
    mc_phase = lguys.Galactocentric[]

    for orbit in df["orbits"]
        @info(orbit)
        icrs = lguys.ICRS(orbit)
        frame = frame_from_params(orbit)

        gc = lguys.transform(lguys.Galactocentric, icrs, frame=frame)
        push!(mc_phase, gc)
    end

    pos = hcat([lguys.position(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0.
    snap = lguys.Snapshot(pos, vel, m)

    return snap#, frames_df_special
end

