using Arya
using CairoMakie
Arya.update_figsize!(390 / 72.27)
Arya.update_fontsize!(12)

CairoMakie.update_theme!(
    linewidth=3,
    markersize=6,
    arrowsize=12,
    Arrows = (;
        linewidth=2,
        arrowsize=9
       ),
    Legend = (;
        patchsize=(18, 6)
   ),
    ErrorScatter = (;
        linewidth=1,
   )
)


function scale_theme_element!(key, scale)
    @info "old $key value = $(theme(key)[])"
    update_theme!(; (; key => scale*theme(key)[])...)
    @info "new $key value = $(theme(key)[])"
end
