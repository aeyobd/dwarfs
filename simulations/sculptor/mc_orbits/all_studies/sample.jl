import LilGuys as lguys
using Distributions

import TOML



obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"

obs_props_all = TOML.parsefile(
    ENV["DWARFS_ROOT"] * "/observations/sculptor/obs_props_all.toml")


obs = lguys.coord_from_file(obs_props_filename)


distributions = Dict{Symbol, Distribution}()
for key in [:ra, :dec, :pmra, :pmdec, :distance, :radial_velocity]
    μs = obs_props_all["$key"]
    σs = obs_props_all["$key" * "_err"]
    distributions[key] = MixtureModel(
        [Normal(μ, σ) for (μ, σ) in zip(μs, σs)]
    )
end

println(distributions)

function sample(N = 100_000)

    obs_special = [obs]

    mc_obs = Vector{lguys.ICRS}(undef, N)

    for i in 1:N
        obs = lguys.ICRS(
            ra = rand(distributions[:ra]),
            dec = rand(distributions[:dec]),
            pmra = rand(distributions[:pmra]),
            pmdec = rand(distributions[:pmdec]),
            distance = rand(distributions[:distance]),
            radial_velocity = rand(distributions[:radial_velocity])
        )
        mc_obs[i] = obs
    end

    mc_obs = vcat(obs_special, mc_obs)

    mc_phase = lguys.transform.(lguys.Galactocentric, mc_obs)
    pos = hcat([lguys.position_of(p) for p in mc_phase]...)
    vel = hcat([lguys.velocity_of(p) for p in mc_phase]...)

    pos ./= lguys.R2KPC
    vel ./= lguys.V2KMS
    vel .*= -1 # reverse velocities to go backwards in time

    m = 0
    snap = lguys.Snapshot(pos, vel, fill(m, N+1))
    return snap
end


if abspath(PROGRAM_FILE) == @__FILE__
    snap = sample()
    lguys.save("initial.hdf5", snap)
end
