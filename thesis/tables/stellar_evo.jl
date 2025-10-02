using LilGuys
using TOML
using OrderedCollections


function calc_props(snap, distance, delta_t=NaN)
    scalars = LilGuys.StellarScalars(snap, delta_t=delta_t)
    df = LilGuys.struct_to_dict(scalars)
    df[:R_h] = LilGuys.half_light_Radius(snap) 
    df[:delta_t_Gyr] = df[:delta_t] * T2GYR
    df[:sigma_v_kms] = df[:sigma_v] * V2KMS
    df[:r_break_arcmin] = LilGuys.kpc2arcmin(df[:r_break], distance)
    return df
end


function print_stellar_props(model_dir, starsname)
    model_dir = joinpath(ENV["DWARFS_ROOT"], "analysis", model_dir)
    stars_dir_in = joinpath(model_dir, "../stars", starsname)

    orbit_props = TOML.parsefile(joinpath(model_dir, "orbital_properties.toml"))
    idx_f = orbit_props["idx_f"]
    df_probs = LilGuys.read_hdf5_table(stars_dir_in * "/probabilities_stars.hdf5")
    probabilities = df_probs.probability

    out = LilGuys.Output(model_dir, weights=probabilities)
    snap_i = out[1]
    snap_f = out[idx_f]
    
    distance = orbit_props["distance_f"]
    delta_t = orbit_props["t_last_peri"] / T2GYR
    d_i = calc_props(snap_i, distance)
    d_f = calc_props(snap_f, distance, delta_t)
    df_out = OrderedDict( "initial" => d_i, "final" => d_f,)

    println("sigma_i", "\t", d_i[:sigma_v_kms])
    println("sigma_f", "\t", d_f[:sigma_v_kms])
    println("frac mass loss", "\t", (d_i[:bound_mass] - d_f[:bound_mass]) / d_i[:bound_mass])
    println("R_h_i", "\t", d_i[:R_h])
    println("R_h_f", "\t", d_f[:R_h])
    println("break radius / kpc", "\t", d_f[:r_break])
    println("break radius / arcmin", "\t", d_f[:r_break_arcmin])
    println()
end


println("sculptor - exp")
print_stellar_props("sculptor/1e7_new_v31_r3.2/orbit_smallperi", "exp2d_rs0.10")
println("sculptor - plummer")
print_stellar_props("sculptor/1e7_new_v31_r3.2/orbit_smallperi", "plummer_rs0.20")

println()

println("sculptor - lmcflyby- exp")
print_stellar_props("sculptor/1e7_new_v25_r2.5/smallperilmc", "exp2d_rs0.11")
println("sculptor - lmcflyby- plummer")
print_stellar_props("sculptor/1e7_new_v25_r2.5/smallperilmc", "plummer_rs0.20")


println()
println("umi - exp")
modeldir = "ursa_minor/1e7_new_v38_r4.0/orbit_smallperi.5"
print_stellar_props(modeldir, "exp2d_rs0.10")
println("umi - plummer")
print_stellar_props(modeldir, "plummer_rs0.20")

