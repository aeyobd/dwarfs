import LilGuys as lguys


include("../sample_utils.jl")


function (@main)(ARGS)
    snap = sample_special_cases()
    lguys.write("initial.hdf5", snap)
end
