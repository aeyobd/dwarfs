import Agama
using LilGuys
using DataFrames
import TOML
import PythonCall


py_agama = Agama.agama[]
q = 0.1


function (@main)(ARGS)
    beta0 = 0.
    r_a = Inf

    if length(ARGS) < 1
        @error "usage $(@__FILE__) starsdirectory"
    end

    df_params = TOML.parsefile("../distribution_function.toml")["DistributionFunction"]
    df_scales = TOML.parsefile("./scales.toml")
    r_scale = df_scales["r_scale"]
    m_scale = df_scales["M_scale"]
    v_scale = df_scales["v_scale"]


    prof = LilGuys.load_profile(ARGS[1] * "/profile.toml")
    if prof isa LilGuys.Exp2D
        dens = Agama.Potential(type="Sersic", sersicIndex=1, scaleRadius = LilGuys.R_h(prof) / r_scale, mass=prof.M / m_scale, axisRatioZ=q)._py
    end

    @info "constructing DFs"
    df_dm = py_agama.DistributionFunction(; LilGuys.dict_to_tuple(df_params)...)
    potential_dm = py_agama.Potential(file="density_iso.ini")

    gridRadius = [0.025, 25]

    grid_kwargs = (; rminSph=gridRadius[1]*0.5, rmaxSph=gridRadius[end], sizeRadialSph=35, lmaxAngularSph=8)

    scm = py_agama.SelfConsistentModel(;
        verbose = false,
        potential = potential_dm,
        components = PythonCall.pylist([py_agama.Component(df=df_dm, density=potential_dm, disklike=false; grid_kwargs...)]),
        grid_kwargs...)

    for i in 1:4
        ϕ0 = scm.potential.potential([0,0,0])
        @info "iteration $i, Φ(0) = $ϕ0"
        scm.iterate()
    end


    df_stars = py_agama.DistributionFunction(type="QuasiSpherical", potential=scm.potential, density=dens, beta0=beta0, r_a=r_a)

    @info "loading snapshot"
    snap = LilGuys.Snapshot("iso_initial.hdf5")

    af = py_agama.ActionFinder(scm.potential)

    posvel_py = Agama.mat2py(vcat(snap.positions / r_scale, snap.velocities / v_scale))

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
