import LilGuys as lguys
import TOML


include("../sample_utils.jl")

obs_props_filename = "observed_properties.toml"

obs_props = TOML.parsefile(obs_props_filename)



function (@main)(ARGS)
    @info("Sampling initial conditions")
    snap = sample(obs_props)

    lguys.write("initial.hdf5", snap)

    @info("Done")
end
