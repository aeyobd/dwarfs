### A Pluto.jl notebook ###
# v0.20.17

using Markdown
using InteractiveUtils

# ╔═╡ 0125bdd2-f9db-11ef-3d22-63d25909a69a
begin
	using Pkg; Pkg.activate()

	FIGDIR = "figures"

	using LilGuys
	using CairoMakie
	using Arya

end

# ╔═╡ f5c22abc-2634-4774-8516-fbd07aa690aa
include("./paper_style.jl")

# ╔═╡ 0e6f0654-dee0-4fee-b9fa-e7fc68de7b91
import FileIO

# ╔═╡ ce3c9118-bcfc-4fa1-8bae-2cafd9052fb8
import TOML

# ╔═╡ dda07928-a282-4a23-ad00-93bb3a8f9e03
module Utils
	include("./utils.jl")
end

# ╔═╡ 5eaf3b50-886e-47ac-9a7c-80d693bc3c17
CairoMakie.activate!(type=:png, px_per_unit=2)

# ╔═╡ 4b115ec9-6848-4d91-a492-f282d1fb1795
lmc = FileIO.load("resources/LMC_4deg_DSS2_color.jpg")

# ╔═╡ 7d9745e3-c413-4bac-b877-21d8ca9adf78
scl = FileIO.load("resources/Scl_1deg_DES-DR2_ColorIRG.jpg")

# ╔═╡ 228c0e35-73bc-4eb7-b006-ed71e0b8d7bd
umi = FileIO.load("resources/umi_2deg_unwise_plus_gaia.png")

# ╔═╡ 67f24990-55fd-412d-aa99-aa69b9c774a2
fornax = FileIO.load("resources/For_1deg_DES-DR2_ColorIRG.jpg")

# ╔═╡ 00c8ab2a-08e8-467e-a9ba-1793d9af8ea5
function plot_scalebar!(width_degrees, bar_degrees, bar_label)
	lines!([0.1, 0.1+bar_degrees / width_degrees], [0.1, 0.1], space=:relative, color=:white, linewidth = theme(:linewidth)[]/2)
	text!(0.1, 0.1, text=bar_label, color=:white, align=(:left, :bottom), space=:relative, fontsize=10, font="TeX Gyre Heros Makie")
	text!
end

# ╔═╡ 2188831f-d8ed-4a37-ad92-bc1f0d90ebbc
function plot_label!(galaxy, fov)
	text!(0.5, 0.9, align=(:center, :center),  text="$galaxy",
		 space=:relative, fontsize=12, font="TeX Gyre Heros Makie", color=:white )
	
	text!(0.9, 0.1, align=(:right, :bottom), text="$fov",
		 space=:relative, fontsize=10, font="TeX Gyre Heros Makie", color=:white)
end

# ╔═╡ 9f38819d-7ec1-4fe0-b3d3-b628b2bf34cb
function plot_ellipse!(galaxy)
	obs_props = TOML.parsefile(joinpath(ENV["DWARFS_ROOT"], "observations", galaxy, "observed_properties.toml"))

	R_h_rel = obs_props["R_h"]

	Utils.ellipse!(R_h_rel, obs_props["ellipticity"], obs_props["position_angle"], color=(:white, 0.5), linewidth=theme(:linewidth)[]/2)
end

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "galaxy_pictures" let
	scalebar = 200
    fig = Figure()

	ax_lmc = Axis(fig[1,1], xreversed=true)
	hidedecorations!()
	fov = 4
	r_arcmin = fov / 2 * 60
	image!((r_arcmin, -r_arcmin), (-r_arcmin, r_arcmin), rotr90(lmc))
	d_lmc = 49.59

	@info LilGuys.arcmin2kpc(0.5*60, d_lmc)
	plot_scalebar!(fov, LilGuys.kpc2arcmin(scalebar / 1e3, d_lmc) / 60, "$scalebar pc")
	plot_label!("LMC", L"1.5\times 10^9\ {L}_\odot")


	ax_boov = Axis(fig[1,2], xreversed=true )
	hidedecorations!()
	d_fornax = 142.56
	@info LilGuys.arcmin2kpc(0.03*60, d_fornax)
	fov = 1
	
	r_arcmin = fov / 2 * 60
	image!((r_arcmin, -r_arcmin), (-r_arcmin, r_arcmin), rotr90(fornax))

	plot_scalebar!(1, LilGuys.kpc2arcmin(scalebar/1e3, d_fornax) / 60, "$scalebar pc")
	plot_label!("Fornax",  L"2\times 10^7\ {L}_\odot")
	plot_ellipse!("fornax")

	

	ax_scl = Axis(fig[2,1], xreversed=true )
	hidedecorations!()
	d_scl = 83.2
	fov = 1
	r_arcmin = fov / 2 * 60
	image!((r_arcmin, -r_arcmin), (-r_arcmin, r_arcmin), rotr90(scl))
	@info LilGuys.arcmin2kpc(0.25*60, d_scl)

	plot_scalebar!(1, LilGuys.kpc2arcmin(scalebar/1e3, d_scl) / 60, "$scalebar pc")
	plot_label!("Sculptor",  L"2\times 10^6\ {L}_\odot")
	plot_ellipse!("sculptor")


	ax_umi = Axis(fig[2,2], xreversed=true )
	hidedecorations!()
		   
	d_umi = 70.1
	fov = 2
	@info LilGuys.arcmin2kpc(0.25*60, d_umi)
		   
	r_arcmin = fov / 2 * 60
	image!((r_arcmin, -r_arcmin), (-r_arcmin, r_arcmin), rotr90(umi))


	plot_scalebar!(2, LilGuys.kpc2arcmin(scalebar/1e3, d_umi) / 60, "$scalebar pc")
	plot_label!("Ursa Minor",  L"3.5\times 10^5\ {L}_\odot")
	plot_ellipse!("ursa_minor")




	rowsize!(fig.layout, 1, Aspect(1, 1.0))

	rowsize!(fig.layout, 2, Aspect(1, 1.0))

	resize_to_layout!(fig)

	rowgap!(fig.layout, 12)
	colgap!(fig.layout, 12)

	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═0e6f0654-dee0-4fee-b9fa-e7fc68de7b91
# ╠═ce3c9118-bcfc-4fa1-8bae-2cafd9052fb8
# ╠═dda07928-a282-4a23-ad00-93bb3a8f9e03
# ╠═f5c22abc-2634-4774-8516-fbd07aa690aa
# ╠═5eaf3b50-886e-47ac-9a7c-80d693bc3c17
# ╠═4b115ec9-6848-4d91-a492-f282d1fb1795
# ╠═7d9745e3-c413-4bac-b877-21d8ca9adf78
# ╠═228c0e35-73bc-4eb7-b006-ed71e0b8d7bd
# ╠═67f24990-55fd-412d-aa99-aa69b9c774a2
# ╠═00c8ab2a-08e8-467e-a9ba-1793d9af8ea5
# ╠═2188831f-d8ed-4a37-ad92-bc1f0d90ebbc
# ╠═9f38819d-7ec1-4fe0-b3d3-b628b2bf34cb
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
