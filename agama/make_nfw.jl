#!/usr/bin/env julia
using PythonCall

using ArgParse
using LilGuys
import TOML

agama = pyimport("agama")
SCRIPT_VERSION = "v0.1.2"


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
        "-b", "--beta"
            help="Velocity anisotropy parameter"
            arg_type=Float64
            default=0.0
        "-a", "--r-a"
            help="Anisotropy scale radius"
            arg_type=Float64
            default=Inf
    end

    args = parse_args(s)


    if args["output"] === nothing
        filename = "halos/nfw_$(args["number"])"

        if args["beta"] != 0.0
            filename = filename * "_beta$(args["beta"])"
        end
        if args["r-a"] != Inf
            filename = filename * "_ra$(args["r-a"])"
        end
        if args["cutoff"] != 100
            filename = filename * "_t$(args["cutoff"])"
        end

        @info "Saving to $filename.hdf5"
        args["output"] = filename * ".hdf5"
    end

    args["number"] = convert(Int, parse(Float64, args["number"]))

    if args["number"] < 1
        throw(ArgumentError("N must be positive"))
    end

    return args
end


function main()
    @info "$(@__FILE__) version $SCRIPT_VERSION"
    @info "LilGuys version $(pkgversion(LilGuys))"

    args = get_args()

    if isfile(args["output"])
        throw(ArgumentError("Output file already exists"))
    end

    @info "args = $args"

    cutoff = args["cutoff"]
    N = args["number"]

    r_s = 1
    M_s = 1
    rho_s = M_s / (4π * r_s^3) 


    pot = agama.Potential(type="Spheroid", densityNorm=rho_s,
        gamma=1, beta=3, alpha=1, scaleRadius=r_s,
        outerCutoffRadius = cutoff*r_s,
        cutoffStrength=1)

    df = agama.DistributionFunction(type="QuasiSpherical", potential=pot,
        beta0=args["beta"], r_a=args["r-a"])

    gm = agama.GalaxyModel(pot, df)

    posvel, mass = gm.sample(N)
    posvel = pyconvert(Matrix{Float64}, posvel)
    posvel = posvel'

    pos = posvel[1:3, :]
    vel = posvel[4:6, :]
    mass = pyconvert(Float64, mass[1])

    snap = Snapshot(pos, vel, mass) 
    LilGuys.write(args["output"], snap)

    halo_kwargs = Dict(
       "profile.TruncNFW" => Dict("trunc" => cutoff, "M_s"=>M_s, "r_s"=>r_s),
       "beta0" => args["beta"],
       "r_a_om" => args["r-a"],
      )

    tomlfilename = splitext(args["output"])[1] * ".toml"

    open(tomlfilename, "w") do f
        TOML.print(f, halo_kwargs)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
