using Arya
using CairoMakie
Arya.update_figsize!(3.54)
Arya.update_fontsize!(10)

CairoMakie.update_theme!(
    linewidth=2,
    markersize=6,
    Legend = (;
        patchsize=(15, 5)
   ),
    ErrorScatter = (;
        linewidth=1,
        markersize=6,
   ),
    Annotation = (;
      style = Ann.Styles.LineArrow(head=Ann.Arrows.Head(length=5)),
      linewidth = 1
    )
)

CairoMakie.activate!(type=:png, px_per_unit=4)

function scale_theme_element!(key, scale)
    @info "old $key value = $(theme(key)[])"
    update_theme!(; (; key => scale*theme(key)[])...)
    @info "new $key value = $(theme(key)[])"
end
