using Arya
using CairoMakie
Arya.update_figsize!(6)
Arya.update_fontsize!(12)

CairoMakie.update_theme!(
    linewidth=3,
    markersize=9,
    Legend = (;
        patchsize=(18, 6)
   ),
    ErrorScatter = (;
        linewidth=1.5,
        markersize=9,
   ),
    Annotation = (;
      style = Ann.Styles.LineArrow(head=Ann.Arrows.Head(length=9)),
      linewidth = 2
    )
)

CairoMakie.activate!(type=:png, px_per_unit=2)

function scale_theme_element!(key, scale)
    @info "old $key value = $(theme(key)[])"
    update_theme!(; (; key => scale*theme(key)[])...)
    @info "new $key value = $(theme(key)[])"
end
