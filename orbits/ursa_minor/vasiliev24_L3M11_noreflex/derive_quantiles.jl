using LilGuys
using PyFITS

using Distributions
using ArgParse
using OrderedCollections


COLUMN_NAMES = ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
time_0 = -9/T2GYR


include("../../orbit_utils.jl")
include("../../derive_quantiles.jl")

function is_bound_to_lmc(df_props, lmc_pot, lmc_centre)
    units = Agama.VASILIEV_UNITS

    is_bound = fill(false, size(df_props, 1))

    for i in eachindex(is_bound)
        pos = [df_props.x_i[i], df_props.y_i[i], df_props.z_i[i]]
        vel = [df_props.v_x_i[i], df_props.v_y_i[i], df_props.v_z_i[i]] / V2KMS
        phi = Agama.potential(lmc_pot, pos, units, t=time_0)

        ke = radii(vel, lmc_centre.velocities[:, 1]) ^ 2 / 2
        is_bound[i] = ke + phi < 0
    end

    return is_bound
end

function is_within_lmc(df_props, lmc_centre, r_max=100)
    units = Agama.VASILIEV_UNITS

    is_bound = fill(false, size(df_props, 1))

    for i in eachindex(is_bound)
        pos = [df_props.x_i[i], df_props.y_i[i], df_props.z_i[i]]

        is_bound[i] = radii(pos, lmc_centre.positions[:, 1]) < r_max
    end

    return is_bound

end

function (@main)(ARGS)
    args = get_args()
    modelname = args["input"]
    df_props = read_fits(joinpath(modelname, "orbital_properties.fits"))


    lmc_pot = Agama.Potential(file="potential_lmc.ini")
    lmc_orbit = get_lmc_orbit(".")
    lmc_cen = LilGuys.resample(lmc_orbit, [time_0])
    println(lmc_cen)


    bound_filter = is_bound_to_lmc(df_props, lmc_pot, lmc_cen)
    distance_filter = is_within_lmc(df_props, lmc_cen, 30)
    @info "fraction bound at end: $(LilGuys.mean(bound_filter))"
    @info "fraction within 30kpc at end: $(LilGuys.mean(distance_filter))"

    df_props = df_props[bound_filter, :]

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
