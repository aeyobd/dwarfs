#!/usr/bin/env julia
using PythonCall

using ArgParse
using LilGuys

agama = pyimport("agama")


function get_args()
    s = ArgParseSettings(
        description="Generate a snapshot of N particles from an NFW profile"
        )

    @add_arg_table s begin
        "output" 
            help="output file. If not provided, will be halos/N_nfw.hdf5 where N is the number of particles"
            default = nothing
        "-N", "--number"
            help="number of particles"
            default="1e4"
        "-t", "--cutoff"
            help="exponential cutoff radius in units of scale radius"
            arg_type=Float64
            default=100
    end

    args = parse_args(s)


    if args["output"] === nothing
        args["output"] = "halos/nfw_$(args["number"]).hdf5"
    end

    args["number"] = convert(Int, parse(Float64, args["number"]))

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

    pos = posvel[1:3, :]
    vel = posvel[4:6, :]
    mass = pyconvert(Float64, mass[1])


    snap = Snapshot(pos, vel, mass) 

    LilGuys.save(args["output"], snap)
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
