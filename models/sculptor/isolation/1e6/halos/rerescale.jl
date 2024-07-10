import LilGuys as lguys

using ArgParse
import TOML



function rescale(snap::lguys.Snapshot, m_scale::Float64, r_scale::Float64)
    v_scale = sqrt(lguys.G*m_scale/r_scale)

    positions = snap.positions * r_scale
    velocities = snap.velocities * v_scale
    masses = snap.masses * m_scale

    return lguys.Snapshot(positions=positions, velocities=velocities, masses=masses, index=snap.index, header=snap.header)
end


function get_args()
    s = ArgParseSettings(description="rescales a snapshot")
    @add_arg_table s begin
        "params"
            help="Parameters file"
            required=true
        "--input"
            help="Input file"
            default="fiducial.hdf5"
        "--initial-params"
            help="Initial parameters file"
            default="fiducial.toml"

        "--output"
            help="Output file"
            default=nothing
    end

    params = ArgParse.parse_args(s)
    if params["output"] === nothing
        params["output"] = splitext(params["params"])[1] * ".hdf5"
    end

    return params
end


function main()
    args = get_args()
    params = TOML.parsefile(args["params"])
    initial_params = TOML.parsefile(args["initial-params"])

    snap = lguys.Snapshot(args["input"])
    prof_old = lguys.NFW(; lguys.dict_to_tuple(initial_params)...)
    prof_new = lguys.NFW(; lguys.dict_to_tuple(params)...)

    m_scale = prof_new.M_s/prof_old.M_s
    r_scale = prof_new.r_s/prof_old.r_s

    snap_rescaled = rescale(snap, m_scale, r_scale)


    println("radius scale: ", r_scale)
    println("mass scale: ", m_scale)

    lguys.save(args["output"], snap_rescaled)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
