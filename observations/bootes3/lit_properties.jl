### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 44b9fa1a-1e3b-11f1-aadc-dbbaf9c22418
begin
	import Pkg;
	Pkg.activate()

	using Arya
	using CairoMakie
end

# ╔═╡ 1cdd4a53-5e25-42c3-9de4-15a5940cb305
using LilGuys

# ╔═╡ 7da5c6fd-1b57-46aa-957c-fb05b5d953e9
±(a, b) = Measurement(a, b)

# ╔═╡ 3a61088e-502c-46ca-b35f-e1ed9de03528
adopted = Dict(
	:radial_velocity => 190 ± 1.8,
	:sigma_v => Measurement(7.5, 1.5, 2.0),
	:pmra => -1.160 ± 0.037,
	:pmdec => -0.88 ± 0.037,
	:ra => 209.5 ± 0.1,
	:dec => 26.7 ± 0.1,
)


# ╔═╡ 1c561ca1-f606-498a-971d-339ed468912d
grillmair09 = Dict(
	:study => "Grillmair2009",
	:ra => 209.3 ± 0,
	:dec => 26.8 ± 0,
	:distance_modulus => 18.35 ± 0.01,
	:pmra => −1.05 ± 0.09,
	:pmdec => −0.95 ± 0.03,
)

# ╔═╡ 80f3cef0-ec82-4558-bc95-2f92469e1848
correnti09 = Dict(
	:study => "Correnti+2009",
	:M_V => -5.8 ± 0.5,
	:ellipticity => 0.5,
	:distance_modulus => 18.58 ± 0.15, 
	:ra => 209.7 ± 1.4/√38,
	:dec => 26.8 ± 0.6/√38,
)

# ╔═╡ 10336f21-2c13-4c00-81eb-86fc4167b344
carlin09 = Dict(
	:study => "Carlin+2009",
	:radial_velocity => 197 ± 3.8,
	:sigma_v => 14 ± 3.2
)


# ╔═╡ 74b809d5-fd34-4b83-badc-ddf92985219c
carlinsand18 = Dict(
	:study => "Carlin+Sand18",
	:radial_velocity => 197.1 ± 3.6,
	:sigma_v => Measurement(10.7, 3.5),
	:distance_modulus => LilGuys.kpc2dm(46.5±2),
	:pmra => −1.14 ± 0.18,
	:pmdec => −0.98 ± 0.20,
)

# ╔═╡ 40c91b00-5bca-47e4-a31b-2f4aef752f9b
massari_helmi18 = Dict(
	:study => "Massari+Helmi2018",
	:pmra => −1.21 ± 0.13,
	:pmdec => −0.92 ± 0.17,
)

# ╔═╡ 498a5787-a2d2-4db4-ae1e-e237fe55e55f
vivas = Dict(
	:study => "vivas+2020",
	:distance_modulus => 18.34 ± 0.19,
)

# ╔═╡ 3bbd9b99-9d6a-482b-9b95-b2e205abb607
LilGuys.dm2kpc(18.34 ± 0.19)

# ╔═╡ 54551965-b8e2-495f-ae7f-1f2cc2661b2c
pace = Dict(
	:study => "pace+2022",
	:pmra => -1.176 ± 0.019,
	:pmdec => -0.890 ± 0.015,
)

# ╔═╡ ff20b3c9-4ceb-40d7-895f-386566750129
geha26 = Dict(
	:study => "Geha+26",
	:radial_velocity => 190 ± 1.8,
	:sigma_v => Measurement(5.2, 1.7, 2.1)
)

# ╔═╡ e234f57c-568d-4097-81ff-16570b1e47c1
myobs = Dict(
	:study => "This work",
	:radial_velocity => 188.9 ± 2.1,
	:sigma_v => Measurement(7.5, 1.5, 2.0)
)

# ╔═╡ dc3d22e8-f7cb-44cd-b0a8-bcd709d962c3
obs = reverse([
	myobs,
	geha26, 
	pace,
	vivas,
	carlinsand18,
	carlin09,
	correnti09,
	grillmair09
])


# ╔═╡ bdc03ed0-ac6e-429e-9d4b-3f1e810a5fe0
function get_properties(obs, key)
	filt = haskey.(obs, key)
	values = [v[key] for v in obs[filt]]
	studies = string.([v[:study] for v in obs[filt]])
	return studies, values
end

# ╔═╡ db9561ea-6cc5-4846-8b0a-58006cd7520b
Arya.value(a::Measurement) = middle.(a)

# ╔═╡ f144cdae-f33a-4a15-81e7-3a31e37aee3f
Arya.err(a::Measurement) = error_interval.(a)

# ╔═╡ 78de8f2b-f62f-4644-b2b7-b19fac4cad8d
function compare_measurements(key, label; units=1, kwargs...)
	fig = Figure()

	x, y = get_properties(obs, key)
	N = length(x)
	y = y / units
	println(y)

	if length(y) == 0
		return fig
	end

	xt = collect(1:N)
	
	ax = Axis(fig[1,1],
		yticks =(xt, x), 
		yminorticksvisible=false,
		# xticklabelrotation=-0π/6,
		ylabel=label, 
		kwargs...
	)

	tight_xticklabel_spacing!(ax)

	errorscatter!(Arya.value.(y), xt, xerror=Arya.err.(y))


	if key ∈ keys(adopted)
		x0 = Arya.value.(adopted[key] / units)
		xl, xh = LilGuys.error_interval.(adopted[key] / units)
		vlines!(x0, color=:black)
		vspan!(x0-xl, x0+xh, color=(:black, 0.3))
	end
	
	fig

end


# ╔═╡ 53952e84-27d5-4ceb-b594-85a5ac94d23a
compare_measurements(:radial_velocity, "radial velocity (km/s)")

# ╔═╡ 15333b02-80df-479b-af9a-9ff5174eae36
compare_measurements(:sigma_v, "velocity dispersion (km/s)")

# ╔═╡ 2801bab1-86f5-4c57-9804-d1c5c230c159
compare_measurements(:distance_modulus, "DM")

# ╔═╡ e1c59c32-8e5d-478b-9191-59d29cd7a42b
compare_measurements(:ra, "RA / deg")

# ╔═╡ a478b91a-3254-4acb-9242-5a48b8db1bea
compare_measurements(:dec, "Dec / deg")

# ╔═╡ b0a6f1af-fa48-4343-81b3-e3d7c6cd3426
compare_measurements(:pmra, "pmra")

# ╔═╡ 718c4805-fa6f-4b61-b1e2-0725cb083b66
compare_measurements(:pmdec, "pmdec")

# ╔═╡ Cell order:
# ╠═44b9fa1a-1e3b-11f1-aadc-dbbaf9c22418
# ╠═1cdd4a53-5e25-42c3-9de4-15a5940cb305
# ╠═7da5c6fd-1b57-46aa-957c-fb05b5d953e9
# ╠═3a61088e-502c-46ca-b35f-e1ed9de03528
# ╠═1c561ca1-f606-498a-971d-339ed468912d
# ╠═80f3cef0-ec82-4558-bc95-2f92469e1848
# ╠═10336f21-2c13-4c00-81eb-86fc4167b344
# ╠═74b809d5-fd34-4b83-badc-ddf92985219c
# ╠═40c91b00-5bca-47e4-a31b-2f4aef752f9b
# ╠═498a5787-a2d2-4db4-ae1e-e237fe55e55f
# ╠═3bbd9b99-9d6a-482b-9b95-b2e205abb607
# ╠═54551965-b8e2-495f-ae7f-1f2cc2661b2c
# ╠═ff20b3c9-4ceb-40d7-895f-386566750129
# ╠═e234f57c-568d-4097-81ff-16570b1e47c1
# ╠═dc3d22e8-f7cb-44cd-b0a8-bcd709d962c3
# ╠═bdc03ed0-ac6e-429e-9d4b-3f1e810a5fe0
# ╠═db9561ea-6cc5-4846-8b0a-58006cd7520b
# ╠═f144cdae-f33a-4a15-81e7-3a31e37aee3f
# ╠═78de8f2b-f62f-4644-b2b7-b19fac4cad8d
# ╠═53952e84-27d5-4ceb-b594-85a5ac94d23a
# ╠═15333b02-80df-479b-af9a-9ff5174eae36
# ╠═2801bab1-86f5-4c57-9804-d1c5c230c159
# ╠═e1c59c32-8e5d-478b-9191-59d29cd7a42b
# ╠═a478b91a-3254-4acb-9242-5a48b8db1bea
# ╠═b0a6f1af-fa48-4343-81b3-e3d7c6cd3426
# ╠═718c4805-fa6f-4b61-b1e2-0725cb083b66
