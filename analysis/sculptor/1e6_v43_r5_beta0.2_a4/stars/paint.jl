import Agama
using LilGuys
using DataFrames
py_agama = Agama.agama[]


function (@main)(ARGS)
    beta0 = 0.2
    r_a = 4.0

    if length(ARGS) < 1
        @error "usage $(@__FILE__) starsdirectory"
    end


    halo = LilGuys.load_profile("../halo.toml")
    pot = LilGuys.AgamaPotential(halo)._py

    prof = LilGuys.load_profile(ARGS[1] * "/profile.toml")
    if prof isa LilGuys.Exp2D
        dens = Agama.Potential(type="Sersic", sersicIndex=1, scaleRadius = LilGuys.R_h(prof), mass=prof.M)._py
    end

    @info "constructing DFs"
    df_stars = py_agama.DistributionFunction(type="QuasiSpherical", potential=pot, density=dens, beta0=beta0, r_a=r_a)
    df_dm = py_agama.DistributionFunction(type="QuasiSpherical", potential=pot, beta0=beta0, r_a=r_a)

    @info "loading snapshot"
    snap = LilGuys.Snapshot("iso_initial.hdf5")

    af = py_agama.ActionFinder(pot)

    posvel_py = Agama.mat2py(vcat(snap.positions, snap.velocities))

    @info "computing actions"
    actions = af(posvel_py)

    @info "computing weights"
    dm_weights = df_dm(actions) |> Agama.py2vec
    stars_weights = df_stars(actions) |> Agama.py2vec


    probs = stars_weights ./ dm_weights

    probs = normalize_probabilities(probs)

    df = DataFrame(
        :index => snap.index,
        :probability => probs,
        :f_s => stars_weights,
        :f_dm => dm_weights,
       )
    sort!(df, :index)

    @info "writing data"
    LilGuys.write_hdf5_table(ARGS[1] * "/probabilities_stars.hdf5", df; overwrite=true)
end


    

function normalize_probabilities(ps)
    N_neg = sum(ps .< 0)
	@info "$N_neg negative probabilities"
    N_nan = sum(isnan.(ps))
    @info "$N_nan NaN probabilities"
	ps[ps .< 0] .= 0
	ps[isnan.(ps)] .= 0

    if sum(ps) == 0
        error("sum of probabilities is zero")
    end

	ps ./= sum(ps)

    return ps
end
