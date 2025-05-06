using CairoMakie
using Arya

CairoMakie.activate!(type=:svg)


Arya.update_figsize!(4)
Arya.update_fontsize!(16)


CairoMakie.update_theme!(
    linewidth=2,
    markersize=6,
    arrowsize=12,
	fonts = (;
		:regular => "TeX Gyre Heros Makie",
		:bold => "TeX Gyre Heros Makie Bold",
		:italic => "TeX Gyre Heros Makie Italic",
		:bold_italic => "TeX Gyre Heros Makie Bold Italic",
	),
    Arrows = (;
        linewidth=2,
        arrowsize=9
    ),
    Legend = (;
        patchsize=(18, 6)
    ),
    ErrorScatter = (;
        linewidth=1,
    ),
)
