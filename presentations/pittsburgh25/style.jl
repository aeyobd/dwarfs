using CairoMakie
using Arya



function scale_theme_element!(key, scale)
    @info "old $key value = $(theme(key)[])"
    update_theme!(; (; key => scale*theme(key)[])...)
    @info "new $key value = $(theme(key)[])"
end

let

    set_theme!(theme_arya(width=860/72, fontsize=40, px_per_unit=1))
    CairoMakie.activate!(type=:svg, px_per_unit=1, pt_per_unit=1)

    scale_theme_element!(:linewidth, 2)
    legend_attr = theme(:Legend)
    legend_attr.margin = legend_attr.padding
    update_theme!(Legend = legend_attr, figure_padding=theme(:fontsize)[])
end
