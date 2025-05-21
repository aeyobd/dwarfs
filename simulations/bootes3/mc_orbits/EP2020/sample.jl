import LilGuys as lguys
import TOML


include("../sample_utils.jl")

obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/bootes3/observed_properties.toml"

obs_props = TOML.parsefile(obs_props_filename)



function (@main)(ARGS)
    @info("Sampling initial conditions")
    snap = sample(obs_props, 100_000)

    lguys.write("initial.hdf5", snap)

    @info("Done")
end
