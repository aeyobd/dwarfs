import TOML
using LilGuys
using CSV, DataFrames
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")


@main function main(args)
    params = read_paramfile("../density_profiles/fiducial.toml")
    pop!(params, "profile_kwargs")

    props = TOML.parsefile("../observed_properties.toml")

    params = LilGuys.dict_to_tuple(params)
    filt_params = GaiaFilterParams(props; params...)

    stars = read_gaia_stars(filt_params)
    members = select_members(stars, filt_params)

    df = members[!, ["source_id", "R_ell"]]

    CSV.write("fornax_members.csv", df)

    return 0
end

