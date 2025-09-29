#!/usr/bin/env julia
using PythonCall

using ArgParse
using LilGuys
import TOML

agama = pyimport("agama")
np = pyimport("numpy")


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
    end

    args = parse_args(s)


    if args["output"] === nothing
        args["output"] = "exp2d_$(args["number"]).hdf5"
    end

    args["number"] = convert(Int, parse(Float64, args["number"]))

    if args["number"] < 1
        throw(ArgumentError("N must be positive"))
    end

    return args
end


function main()
    args = get_args()

    N = args["number"]

    halo = LilGuys.Exp2D(M=1.0, R_s=1.0)

    function ρ(x)
        xj = pyconvert(Matrix{Float64}, x)'
        r = radii(xj)
        rho =  np.array(LilGuys.density.(halo, r))
        return rho
    end

    pot = agama.Potential(type="Multipole", 
        density = pyfunc(ρ),
        symmetry="s",
        gridSizeR=200, rmin=0.001, rmax=30
       )

    df = agama.DistributionFunction(type="QuasiSpherical", potential=pot)
    gm = agama.GalaxyModel(pot, df)

    posvel, mass = gm.sample(N)
    posvel = pyconvert(Matrix{Float64}, posvel)
    posvel = posvel'

    pos = posvel[1:3, :]
    vel = posvel[4:6, :]
    mass = pyconvert(Float64, mass[1])


    snap = Snapshot(pos, vel, mass) 

    LilGuys.write(args["output"], snap)

    halo_kwargs = Dict( "Exp2D" => LilGuys.struct_to_dict(halo))

    tomlfilename = args["output"][1:end-5] * ".toml"
    open(tomlfilename, "w") do f
        TOML.print(f, halo_kwargs)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
