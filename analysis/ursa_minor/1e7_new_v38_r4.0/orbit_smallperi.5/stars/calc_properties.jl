using LilGuys
using TOML
using OrderedCollections


model_dir = ".."
starsname = ARGS[1]

stars_dir_in = joinpath(model_dir, "../stars", starsname)

@info "loading info"
orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))
idx_f = orbit_props["idx_f"]
df_probs = LilGuys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")
probabilities = df_probs.probability



@info "loading snapshots"
out = LilGuys.Output(model_dir, weights=probabilities)
snap_i = out[1]
snap_f = out[idx_f]


function calc_props(snap, delta_t=NaN)
    scalars = LilGuys.StellarScalars(snap, delta_t=delta_t)
    df = LilGuys.struct_to_dict(scalars)
    df[:R_h] = LilGuys.half_light_Radius(snap) 
    df[:delta_t_Gyr] = df[:delta_t] * T2GYR
    df[:sigma_v_kms] = df[:sigma_v] * V2KMS
    df[:r_break_arcmin] = LilGuys.kpc2arcmin(df[:r_break], orbit_props["distance_f"])
    return df
end


outname = joinpath(starsname, "properties.toml")
open(outname, "w") do f
    d_i = calc_props(snap_i)
    d_f = calc_props(snap_f, orbit_props["t_last_peri"] / T2GYR)
    df_out = OrderedDict( "initial" => d_i, "final" => d_f,)


    TOML.print(f, df_out)
    TOML.print(df_out)
end
