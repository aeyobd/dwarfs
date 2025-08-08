using LilGuys
using PyFITS

using Distributions
using ArgParse
using OrderedCollections


COLUMN_NAMES = ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]


include("orbit_utils.jl")

function get_args()
    s = ArgParseSettings(
        description="""Calculates the orbits
        returing essential information.
""",
    )

    @add_arg_table s begin
        "input"
            help="input directory, should contain agama_potential.ini"
        "output"
            help="output file. defaults to median_properties.toml"
        "--galaxy"
            help = "the name of the galaxy"
        "--key"
            help = "the key in orbital_properties to take quantile filters of"
            default = "pericentre"
        "--p-value"
            help = "probability value for quantiles"
            default = cdf(Normal(), -3)
            arg_type = Float64
    end

    args = parse_args(s)


    if isnothing(args["output"])
        args["output"] = joinpath(dirname(args["input"] * "/"), "median_properties.toml")
    end

    return args
end



function median_residual(df_props, obs_props)
    df = OrderedDict()
	for key in COLUMN_NAMES
        x = obs_props[key]
        err = LilGuys.get_uncertainty(obs_props, key)
		xs = df_props[!, key]
		md = median(xs)
		res = (md - x ) / err
        df["$(key)_relative_residual"] = res
	end

    return df
end


function median_properties(df_props)
    df = OrderedDict()
	for key in COLUMN_NAMES
		x = df_props[!, key]
		md = median(x)
		err = std(x) / sqrt(length(x))
        df["$(key)_median"] = md
        df["$(key)_sem"] = err
	end

    return df
end


function quantile_values(df, obs_props)
    return merge(
         median_properties(df),
         median_residual(df, obs_props),
    )
end



function (@main)(ARGS)
    args = get_args()
    modelname = args["input"]
    df_props = read_fits(joinpath(modelname, "orbital_properties.fits"))
    obs_props = get_obs_props(args["input"], args["galaxy"])
    p_value = args["p-value"]
    quantiles = [0.5, p_value, 1-p_value]

    orbit_labels = ["mean", "smallperi", "largeperi"]

    key = args["key"]

    peris = df_props[!, key]

    median_df = OrderedDict()

    for (label, q) in zip(orbit_labels, quantiles)
        if q == 0.5
            df = df_props
        elseif q < 0.5
            cut = quantile(peris, 2q)
            df = df_props[peris .< cut, :]
        elseif q > 0.5
            cut = quantile(peris, 2q - 1)
            df = df_props[peris .> cut, :]
        end

        medians = quantile_values(df, obs_props)
        medians[key] = median(df[!, key])
        medians["quantile"] = q
        medians["n_samples"] = size(df, 1)


        if !isapprox(median(df[!, key]), quantile(peris, q), rtol=1e-2) 
            @warn "internal consistency with quantiles failed, maybe small sample size"
        end

        median_df[label] = medians
    end

    peri_qs = LilGuys.quantile(peris, quantiles)

    open(args["output"], "w") do f
        TOML.print(f, median_df)
    end
end
