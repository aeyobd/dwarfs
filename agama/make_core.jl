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
        "-t", "--cutoff"
            help="exponential cutoff radius in units of scale radius"
            arg_type=Float64
            default=100
        "-c", "--core"
            help="core radius in units of scale radius"
            arg_type=Float64
            default=0
    end

    args = parse_args(s)


    if args["output"] === nothing
        args["output"] = "halos_cored/cnfw_$(args["number"])_rc_$(args["core"]).hdf5"
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
    r_c = args["core"]

    scaleRadius = 1
    mass = 1
    rho0 = mass / (4π * scaleRadius^3) 

    halo = LilGuys.CoredNFW(M_s=1, r_s=1, r_c=r_c, r_t=cutoff)

    function ρ(x)
        xj = pyconvert(Matrix{Float64}, x)'
        r = calc_r(xj)
        rho =  np.array(calc_ρ.(halo, r))
        return rho
    end

    pot = agama.Potential(type="Multipole", 
        density = pyfunc(ρ),
        symmetry="s",
        gridSizeR=200, rmin=0.001, rmax=5*cutoff
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

    LilGuys.save(args["output"], snap)

    halo_kwargs = Dict(
       "CoredNFW" => Dict("r_c" => r_c, "r_t" => cutoff, "M_s"=>1, "r_s"=>1)
      )

    tomlfilename = args["output"][1:end-5] * ".toml"
    open(tomlfilename, "w") do f
        TOML.print(f, halo_kwargs)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
