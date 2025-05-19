import LilGuys as lguys
import TOML
using PyFITS

include("../sample_utils.jl")


obs_props_filename = ENV["DWARFS_ROOT"] * "/observations/sculptor/observed_properties.toml"
obs_props = TOML.parsefile(obs_props_filename)


function (@main)(ARGS)
    snap, frames_df = sample(obs_props, frame_distributions)
    write_fits("gc_frames.fits", frames_df)
    lguys.write("initial.hdf5", snap)
end
