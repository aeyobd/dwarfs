import LilGuys as lguys
using Distributions

import TOML

using DataFrames

obs = lguys.ICRS(
    ra = 15.0385, # M 2012 
    dec = -33.7056,
    pmra = 0.085,
    pmdec = -0.133,
    distance = 86,
    radial_velocity = 111.6,

   )

err = lguys.ICRS(
    ra = 0.,
    dec = 0.,
    pmra = 0.028,
    pmdec = 0.028,
    distance = 6,
    radial_velocity = 0.2,
   )


frame_distributions = Dict{Symbol, Distribution}(
    :d => Normal(8.0, 0.5),
    :ra => Normal(266.4051, 0.0003),
    :dec => Normal(-28.936175, 0.000003),
    :z_sun => Normal(0.0, 0.0),
    :v_sun_x => Normal(11.1, 1.3),
    :v_sun_y => Normal(230.2, 6),
    :v_sun_z => Normal(7.25, 0.6),
   )

mean_frame = lguys.GalactocentricFrame(
    d = 8.0,
    ra = 266.4051,
    dec = -28.936175,
    z_sun = 0.0,
    v_sun = [11.1, 230.2, 7.25]
   )

function sample(N = 100_000)

    obs_special = [obs,
        lguys.ICRS(ra=obs.ra, dec=obs.dec, distance=82.3,
                   pmra=0.125, pmdec=-0.164, radial_velocity=111.8),
                  ]

    frames_special = [mean_frame,
        lguys.GalactocentricFrame(d=8.2, ra=266.4051, dec=-28.936175, z_sun=0.0,
                                v_sun=[10.9, 213.3, 7.9]),
                      ]

    mc_obs = [lguys.rand_coord(obs, err) for _ in 1:N]
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

    # add in special orbits
    mc_obs = vcat(obs_special, mc_obs)

    frames = vcat(frames_special, frames)

    frames_df_special = DataFrame(Dict(
        :d => [frame.d for frame in frames_special],
        :ra => [frame.ra for frame in frames_special],
        :dec => [frame.dec for frame in frames_special],
        :z_sun => [frame.z_sun for frame in frames_special],
        :v_sun_x => [frame.v_sun[1] for frame in frames_special],
        :v_sun_y => [frame.v_sun[2] for frame in frames_special],
        :v_sun_z => [frame.v_sun[3] for frame in frames_special]
    ))

    frames_df = vcat(frames_df_special, frames_df)

    # transform to phase space coordinates
    mc_phase = [lguys.transform(lguys.Galactocentric, o, frame=frame) for (o, frame) in zip(mc_obs, frames)]
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))

    return snap, frames_df
end


function (@main)(ARGS)
    snap, frames_df = sample()
    lguys.write_fits("gc_frames.fits", frames_df)
    lguys.save("initial.hdf5", snap)
end
