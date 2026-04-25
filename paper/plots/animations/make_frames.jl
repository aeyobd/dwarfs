include("../picture_utils.jl")
using LilGuys
using HDF5
import TOML


function (@main)(ARGS)
    idx_i = parse(Int, get(ARGS, 1, "1"))
    idx_f = parse(Int, get(ARGS, 2, "-1"))


    time_f = TOML.parsefile(joinpath("model", "orbital_properties.toml"))["t_f_gyr"] / T2GYR

    stars = LilGuys.read_hdf5_table(joinpath("stars", "probabilities_stars.hdf5"))
    out = Output("model", weights=stars.probability)

    snap_i = out[1]
    colorrange_dm = (-5, 0) .+ log10(maximum(get_zoom_histogram(snap_i)[2]))
    colorrange_stars = (-5, 0) .+ log10(maximum(get_zoom_histogram(snap_i, weights=snap_i.weights)[2]))

    if idx_f == -1
        idx_f = length(out)
    end
    for i in eachindex(out)[idx_i:idx_f]
        @info "animating frame $i"
        fig = plot_frame(out[i], colorrange_dm=colorrange_dm, colorrange_stars=colorrange_stars)
        time = out.times[i] - time_f

        add_time!(fig.content[1], time * T2GYR)
        Makie.save(joinpath("animation", "frame_$(i).png"), fig)
    end

end
