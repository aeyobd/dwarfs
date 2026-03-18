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

# ╔═╡ 10336f21-2c13-4c00-81eb-86fc4167b344
carlin09 = Dict(
	:study => "Carlin+2009",
	:radial_velocity => 197 ± 3.8,
	:sigma_v => 14 ± 3.2
)


# ╔═╡ e234f57c-568d-4097-81ff-16570b1e47c1
myobs = Dict(
	:study => "This work",
	:radial_velocity => 188.9 ± 2.1,
	:sigma_v => Measurement(7.5, 1.5, 2.0)
)

# ╔═╡ ff20b3c9-4ceb-40d7-895f-386566750129
geha26 = Dict(
	:study => "Geha+26",
	:radial_velocity => 190 ± 1.8,
	:sigma_v => Measurement(5.2, 1.7, 2.1)
)

# ╔═╡ 74b809d5-fd34-4b83-badc-ddf92985219c
carlinsand18 = Dict(
	:study => "Carlin+Sand18",
	:radial_velocity => 197.1 ± 3.6,
	:sigma_v => Measurement(10.7, 3.5)
)

# ╔═╡ dc3d22e8-f7cb-44cd-b0a8-bcd709d962c3
obs = reverse([
	myobs,
	geha26, 
	carlinsand18,
	carlin09
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


	# if true #key ∈ keys(adopted)
	# 	x0 = Arya.value.(adopted[key] / units)
	# 	xe = Arya.err.(adopted[key] / units)
	# 	hlines!(x0, color=:black)
	# 	hspan!(x0-xe, x0+xe, color=(:black, 0.3))
	# end
	
	fig

end


# ╔═╡ 53952e84-27d5-4ceb-b594-85a5ac94d23a
compare_measurements(:radial_velocity, "radial velocity (km/s)")

# ╔═╡ 15333b02-80df-479b-af9a-9ff5174eae36
compare_measurements(:sigma_v, "velocity dispersion (km/s)")

# ╔═╡ Cell order:
# ╠═44b9fa1a-1e3b-11f1-aadc-dbbaf9c22418
# ╠═1cdd4a53-5e25-42c3-9de4-15a5940cb305
# ╠═7da5c6fd-1b57-46aa-957c-fb05b5d953e9
# ╠═10336f21-2c13-4c00-81eb-86fc4167b344
# ╠═e234f57c-568d-4097-81ff-16570b1e47c1
# ╠═ff20b3c9-4ceb-40d7-895f-386566750129
# ╠═74b809d5-fd34-4b83-badc-ddf92985219c
# ╠═dc3d22e8-f7cb-44cd-b0a8-bcd709d962c3
# ╠═bdc03ed0-ac6e-429e-9d4b-3f1e810a5fe0
# ╠═db9561ea-6cc5-4846-8b0a-58006cd7520b
# ╠═f144cdae-f33a-4a15-81e7-3a31e37aee3f
# ╠═78de8f2b-f62f-4644-b2b7-b19fac4cad8d
# ╠═53952e84-27d5-4ceb-b594-85a5ac94d23a
# ╠═15333b02-80df-479b-af9a-9ff5174eae36
