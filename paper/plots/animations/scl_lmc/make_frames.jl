using LilGuys
using HDF5
import TOML

include("../../picture_utils.jl")
include(joinpath(ENV["DWARFS_ROOT"], "orbits/orbit_utils.jl"))


function (@main)(ARGS)
    idx_i = parse(Int, get(ARGS, 1, "1"))
    idx_f = parse(Int, get(ARGS, 2, "-1"))

    bins = LinRange(-300, 300, 1001)

    lmc_orbit = get_lmc_orbit(joinpath(ENV["DWARFS_ROOT"], "orbits/sculptor/vasiliev24_L3M11")) 

    time_f = TOML.parsefile(joinpath("model", "orbital_properties.toml"))["t_f_gyr"] / T2GYR

    stars = LilGuys.read_hdf5_table(joinpath("stars", "probabilities_stars.hdf5"))
    out = Output("model", weights=stars.probability)
    lmc_orbit = LilGuys.resample(lmc_orbit, out.times)

    snap_i = out[1]
    colorrange_dm = (-5, 0) .+ log10(maximum(get_zoom_histogram(snap_i)[2]))
    colorrange_stars = (-5, 0) .+ log10(maximum(get_zoom_histogram(snap_i, weights=snap_i.weights)[2]))


    if idx_f == -1
        idx_f = length(out)
    end
    for i in eachindex(out)[idx_i:idx_f]
        @info "animating frame $i"
        fig = plot_frame(out[i], colorrange_dm=colorrange_dm, colorrange_stars=colorrange_stars, bins=bins)
		scatter!(lmc_orbit.positions[2, i], lmc_orbit.positions[3, i], color=COLORS[3])

        time = out.times[i] - time_f

        add_time!(fig.content[1], time * T2GYR)
        Makie.save(joinpath("animation", "frame_$(i).png"), fig)
    end

end
