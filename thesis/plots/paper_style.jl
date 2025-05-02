using Arya
using CairoMakie
Arya.update_figsize!(390 / 72.27)
Arya.update_fontsize!(12)

CairoMakie.update_theme!(
    linewidth=2,
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
