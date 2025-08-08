using Arya
using Makie
import LilGuys as lguys


ms_special = @lift 2*$(theme(:markersize))

coord_labels = Dict(
	:ra => "ra / degrees",
	:dec => "dec / degrees",
	:pmra => L"$\mu_{\alpha *}$ / mas\,yr$^{-1}$",
	:pmdec => L"$\mu_{\delta}$ / mas\,yr$^{-1}$",
	:radial_velocity => L"$v_\textrm{los}$ / km\,s$^{-1}$",
	:distance => "distance / kpc",

    :pericentre => "pericentre / kpc",
    :apocentre => "apocentre / kpc",
    :pericentre_lmc => "peri-LMC / kpc",
    :apocentre_lmc => "apo-LMC / kpc",
    :time_last_peri => "time of last peri",
    :time_last_peri_lmc => "time last LMC peri",
)



function plot_param_hist(df_props, special_props, y_key)
	fig = Figure()
	ax = Axis(fig[1, 1],
        xlabel = coord_labels[y_key],
		ylabel = "count"
	)

    bins, counts, err = lguys.histogram(df_props[!, y_key])
	scatter!(lguys.midpoints(bins), counts)


    for i in 1:size(special_props, 1)
        scatter!(special_props[i, y_key], 0, color=COLORS[i], label=special_props.label[i], markersize=ms_special, marker=:star5)
	end
	
	axislegend()
	fig
end


function plot_correlations(df_props, special_props, y_key=:pericentre;
        Nmax = 10_000,
        ax_idx = Dict(
            :pmra => [1, 1],
            :pmdec => [1, 2],
            :distance => [2, 1],
            :radial_velocity => [2, 2],
            :ra => [3,1],
            :dec => [3,2],
        )
    )

    Nmax = min(Nmax, size(df_props, 1))

    fig = Figure()

    ax_kwargs = Dict(
        :xgridvisible => false,
        :ygridvisible => false,
        :ylabel => coord_labels[y_key], 
        :width => 2*72,
        :height => 2*72,
    )

    plot_kwargs = Dict(
        :color => :black,
        :alpha => 0.1,
        :markersize => 1,
    )

    orbit_points_kwargs = Dict(
        :alpha => 1,
        :markersize => 10,
    )

    axes = Dict()

    for (sym, idx) in ax_idx
        ax = Axis(fig[idx...];
			xlabel=coord_labels[sym],
			ax_kwargs...
		)

	    x = df_props[1:Nmax, sym]
		y = df_props[1:Nmax, y_key]
		scatter!(x, y; plot_kwargs...)
		
        for (i, row) in enumerate(eachrow(special_props))
            x = row[sym]
            y = row[y_key]
            label = row["label"]
            scatter!(x, y; label=label, color=COLORS[i], orbit_points_kwargs...)
		end

        axes[sym] = ax
	end


	linkyaxes!(fig.content...)

    axislegend(axes[:ra])


    resize_to_layout!(fig)
	fig
end
