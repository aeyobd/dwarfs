using PythonCall

using ArgParse
using LilGuys

agama = pyimport("agama")


function get_args()
    s = ArgParseSettings()

    @add_arg_table s begin
        "output" 
            help="output file"
            default = "nfw.hdf5"
        "-v", "--verbose"
            help="verbose output"
            action="store_true"
        "-N", "--number"
            help="number of particles"
            arg_type=Float64
            default=1e4
        "-t", "--cutoff"
            help="cutoff radius"
            arg_type=Float64
            default=100
    end

    args = parse_args(s)
    args["number"] = convert(Int, args["number"])

    if args["number"] < 1
        throw(ArgumentError("N must be positive"))
    end

    return args
end


function main()
    args = get_args()

    cutoff = args["cutoff"]
    N = args["number"]

    scaleRadius = 1
    mass = 1
    rho0 = mass / (4π * scaleRadius^3) 


    pot = agama.Potential(type="Spheroid", densityNorm=rho0,
        gamma=1, beta=3, alpha=1, scaleRadius=scaleRadius,
        outerCutoffRadius = cutoff*scaleRadius,
        cutoffStrength=1)

    df = agama.DistributionFunction(type="QuasiSpherical", potential=pot)
    gm = agama.GalaxyModel(pot, df)

    posvel, mass = gm.sample(N)
    posvel = pyconvert(Matrix{Float64}, posvel)
    posvel = posvel'
    println(size(posvel))

    pos = posvel[1:3, :]
    vel = posvel[4:6, :]
    mass = pyconvert(Float64, mass[1])


    snap = Snapshot(pos, vel, mass) 

    LilGuys.save(args["output"], snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
