using ArgParse
using LilGuys
import TOML
include(ENV["DWARFS_ROOT"] * "/utils/gaia_filters.jl")


function get_args()
    s = ArgParseSettings(description="applies the given filters to a data file and computes the stellar profile")

    @add_arg_table s begin
        "input"
            help = "path to TOML file for filter & profile settings"
            required = true
        "--quantile", "-q"
            arg_type = Float64
            default = 0.5
            help = "quantile to use"
    end

    args = parse_args(s)

    return args
end

function main()
    args = get_args()

    params = read_paramfile(args["input"])

    profile_kwargs = pop!(params, "profile_kwargs", Dict())

    props = TOML.parsefile(dirname(params["filename"]) * "/../observed_properties.toml")

    params = LilGuys.dict_to_tuple(params)
    filt_params = GaiaFilterParams(props; params...)

    stars = read_gaia_stars(filt_params)
    members = select_members(stars, filt_params)

    q = args["quantile"]
    q_val = LilGuys.quantile(members.G, q)
    println("$q quantile = $q_val")
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

