import LilGuys as lguys
import TOML


include("../sample_utils.jl")

obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/bootes3/observed_properties.toml"

obs_props = TOML.parsefile(obs_props_filename)
obs_props["distance_modulus_em"] = 0.095
obs_props["distance_modulus_ep"] = 0.095
# corresponds to a 2kpc distance error



function (@main)(ARGS)
    @info("Sampling initial conditions")
    snap = sample(obs_props)

    lguys.write("initial.hdf5", snap)

    @info("Done")
end
