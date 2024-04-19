### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ db8a38d4-fae2-11ee-18e3-a189bdef1759
begin
	import Pkg; Pkg.activate()

	import LilGuys as lguys

	using Plots
	using Arya

	using CSV, DataFrames
end

# ╔═╡ dc7c775a-606b-4fb4-a688-4c858c31ca20
using FFTW

# ╔═╡ 04b79ff4-1c47-44bc-a7d1-e3e119dccd99
begin 
	names = ["fiducial", "light", "heavy", "flat", "steep"]

	orbits = Dict(name=> CSV.read("$name/nbody_orbit.csv", DataFrame) for name in names)

	densities = Dict(name=> CSV.read("$name/density_prof.csv", DataFrame) for name in names)
end

# ╔═╡ 681dc280-e68a-4e52-990b-88089e1803e4
function plot_orbit!(orbit; kwargs...)
	R = @. sqrt(orbit[:, "xs"]^2 + orbit[:, "ys"]^2)
	plot!(R, orbit[:, "zs"]; kwargs...)
end

# ╔═╡ 1bdea423-3869-4746-9655-ac7c43aaf652
densities["fiducial"]

# ╔═╡ a499624d-3e0e-445d-a4c8-69be93205d42
p = Ref{Plots.Plot}()

# ╔═╡ d6166e1d-501b-4744-b392-e389b02d7640
begin 
	p[] = plot(xlabel=L"\log r\,/\,\rm pc", ylabel=L"\log \rho\,/\,\rm M_\odot\,pc^{-3}", size=(335, 250))
	for name in names
		plot!(log10.(densities[name].r * 1e3), log10.(densities[name].rho * lguys.M0 / 1e9), label="")
	end
	p[]
end

# ╔═╡ 15ff24b6-533d-454d-acf8-b3e0f25fed39
function plot_obs2!(obs; kwargs...)
	ra = [o.ra for o in obs]
	dec = [o.dec for o in obs]
	ra[ra .< 0] .+= 360
	plot!(ra, dec; kwargs...)
end

# ╔═╡ e7429c5e-216b-4bf6-8dd7-61754af54523
function plot_obs!(obs; kwargs...)
	ra = [o.ra for o in obs]
	dec = [o.dec for o in obs]
	ra[ra .< 100] .+= 360
	scatter!(ra, dec; kwargs...)
end

# ╔═╡ 8c62cdef-4b48-4d96-add6-05fc6e370372
jax = CSV.read("../jax_membs.tsv", DataFrame)

# ╔═╡ 8104512d-c99c-433a-8afc-db949740f3ad
Arya.COLORS[9]

# ╔═╡ 0d2c4538-16e7-41d3-b8d8-bd4872f15026
begin 
	p[] = plot(
		xlims=(-1.5, 0.1),
		ylims=(0, 65),
		xlabel="time / Gyr",
		ylabel="\$r\$ / kpc",
		legend_position=:bottomleft
	)

	for name in names
		plot!(orbits[name].times * lguys.T0, orbits[name].rs, label=name)
	end
	savefig("r_t.pdf")
	p[]
end

# ╔═╡ f994c87e-15df-421d-9282-361d1bb55f16
pythonplot()

# ╔═╡ 124e0435-ccad-4918-bebc-762fa846ee69
Arya.set_mpl_theme()

# ╔═╡ dbc15de6-8c9b-4e25-92dc-c02a68976355
dpi=500

# ╔═╡ 4e103842-f084-4a53-a5ae-c149447744a1
Arya.half_width

# ╔═╡ 1fadb10f-b160-41f7-99b3-25e654d0a674
begin 
	p[] = plot(
		ylabel="bound fraction",
		xlabel="time / Gyr",
		size= 100 .* (3.35, 2.5),
		dpi=dpi,
		legend_position=:bottomleft
	)

	for name in names
		plot!(orbits[name].times * lguys.T0, orbits[name].f, label=name)
	end

	savefig("fraction_versus_time.pdf")
	p[]
end

# ╔═╡ ed11b804-49f8-463b-84fd-7b3c18a0c5f6
orbits["heavy"].times * lguys.T0

# ╔═╡ 3813d261-6890-44aa-9560-c2758745eec4
idx_0s = Dict(name=>findfirst(orbits[name].times .> 0) for name in names)

# ╔═╡ f9b3062d-cae6-4eaf-b9cd-6dc5a79f14ce
begin 
	obsp = Dict{String,Any}()
		
	for name in names
		i_0 = idx_0s[name]
		o = orbits[name][i_0-100:i_0, :]
		
		particles_gc = [lguys.Galactocentric(o[i, "xs"], o[i, "ys"], o[i, "zs"],  0,0,0) for i in 1:length(o.xs)]
	    obs = lguys.transform.(lguys.Observation, particles_gc)
		obsp[name] = obs
	end
end

# ╔═╡ 6ae8ec6c-dc20-497e-9fd2-a8e4cffd920a
begin 
	p[] = plot()
	
	for name in names
		plot_obs2!(obsp[name], label=name, ms=1)
	end
	
	p[]
end

# ╔═╡ 97d8f567-bcbe-4d24-a609-3f2002b71b29
begin
	p[] = plot(aspect_ratio=1)

	for name in names
		i0 = idx_0s[name]
		plot_orbit!(orbits[name][i0-100:i0, :], label=name)
	end
	p[]
end

# ╔═╡ b90632da-680c-49d5-9903-b00191a4f03b
begin 
	p[] = plot(aspect_ratio=1
	)

	for name in names
		i0 = idx_0s[name]
		o = orbits[name][i0-200:i0, :]
		plot!(o.xs , o.ys, label=name)
	end
	
	p[]
end

# ╔═╡ b39e1baf-2de0-4539-9957-200679fcd38d
idx_0s

# ╔═╡ a3acf540-ac34-49f3-8594-6d8c4718c83e
snaps = Dict(name=>lguys.Snapshot("$name/out/snapshot_$(idx_0s[name] -2).hdf5") for name in names)

# ╔═╡ ac0b3eac-015d-4a96-919e-fc90ddb575bf
begin 
	obss = Dict{String,Any}()
		
	for name in names
		snap = snaps[name]
		
		particles_gc = [lguys.Galactocentric(snap.positions[:, i]..., (lguys.V0 * snap.velocities[:, i])...) for i in 1:length(snap)]
	    obs = lguys.transform.(lguys.Observation, particles_gc)
		obss[name] = obs
	end
end

# ╔═╡ 8ff01570-4508-49db-869d-68fab24b3fc0
begin
	p[] = plot(xlabel=L"\rm RA\, / \, ^{\circ}", ylabel=L"\rm Dec\, / \, ^{\circ}", xlims=(150, 250), ylims=(10, 40))
	scatter!(jax.RA, jax.DEC, label="Jensen+2021", ms=2)
	
	plot_obs!(obss["fiducial"], label="", ms=1, alpha=0.3)

	savefig("stream_obs.pdf")
	p[]
end

# ╔═╡ 5136a85e-5f19-4604-b1d1-5ef7e04bdb29
begin 
	p[] = plot(xlabel=L"\rm RA\, / \, ^{\circ}", ylabel=L"\rm Dec\, / \, ^{\circ}" )
	
	for (i, name) in enumerate(names)
		plot_obs!(obss[name], label=name, ms=1, color=Arya.COLORS[i])
	end

	scatter!(jax.RA, jax.DEC, label="Jensen+2021", ms=2, cmap=Arya.COLORS[9])

	savefig("stream_allsky.pdf")
	p[]
end

# ╔═╡ 38d1adb5-a767-429f-bf6a-644e23bb0655
snaps_i = Dict(name=>lguys.Snapshot("$name/out/snapshot_001.hdf5") for name in names)

# ╔═╡ fcbbfb05-f761-42c2-a6e6-1121687b0a62
function plot_snap(snap)
	scatter(snap.positions[1, :], snap.positions[2, :], alpha=0.1, ms=1, label="")
end

# ╔═╡ 874346e2-d285-42b7-aa94-c34d139cb2ee
gr()

# ╔═╡ db5ad747-2e70-4861-a260-7bf0496eefd9
plot_snap(snaps_i["light"])

# ╔═╡ c30590cb-1d28-438d-8c98-31e752303fee
plot_snap(snaps_i["fiducial"])

# ╔═╡ 4d71c42f-7bfe-436b-aae0-9ffe75cbfc81


# ╔═╡ 8ae6acb2-c4fa-4718-9b44-1937d04be1f0
plot_snap(snaps_i["heavy"])

# ╔═╡ 368483d0-84d4-48ef-8f61-95ac001a587f
single_size=100 .* (3.35, 2.5)

# ╔═╡ de6d5e1e-2937-4199-9274-a0bc49fb8d8a
theme(:arya;Arya.make_font_settings(12)..., size=single_size, dpi=300)

# ╔═╡ 6807c973-c193-4f66-9a1e-a569f3ad872f
begin 
	p[] = plot_snap(snaps["fiducial"],)
	plot!(xlabel=L"x\,/\, \rm kpc", ylabel=L"y\,/\, \rm kpc",  size=single_size)
	savefig("xy_scatter.pdf")
	p[]
end

# ╔═╡ 6fecdb3a-9945-4027-ba21-07bd6fec9a36
plot_snap(snaps["light"])

# ╔═╡ c649a8a6-2ede-4988-ba69-c6df360931a5
plot_snap(snaps["heavy"])

# ╔═╡ cf2f33a0-d940-4ddb-bc8c-592af2dff0a7
plot_snap(snaps["steep"])

# ╔═╡ 6de087be-fa5a-4f2e-ab5d-189385dc2676
plot_snap(snaps["flat"])

# ╔═╡ 41b5f8c4-f514-4c83-a46e-2c4e7fc8b492
o = orbits["fiducial"]

# ╔═╡ cc81994f-3944-4e01-bb8e-cfc2982c800c
begin
	plot()

	plot!(o.times, o.f_i, label="inner")
	plot!(o.times, o.f_m, label="middle")
	plot!(o.times, o.f_o, label="outer")
end

# ╔═╡ 3d02ec0a-7376-4530-a66f-7297454d98e1
begin
	plot()
	plot!(o.times, o.M_3r ./ o.M_3r[1])
end

# ╔═╡ d36a3b1c-d224-4960-9fa6-aab119b8ab61
begin 

end

# ╔═╡ f201f6c4-f2dc-4e28-91b7-cd3f5ce871b2
function compute_period(times, radii)
	freq =  1/lguys.T0*1/lguys.mean(diff(times)) 
	N = length(times)
	freqs =  rfftfreq(N, freq)
	periods = 1 ./ freqs
	
	f = abs.(rfft(radii))
	idx = argmax(f[2:end]) + 1 # ignore zero point
	return periods[idx]
end

# ╔═╡ 05bd3bd4-93dc-48d8-ac97-ecddc9390fca
for name in names
	println(name)
	r = orbits[name].rs
	println("peri\t ", minimum(r))

	println("apo \t ", maximum(r))
	T = compute_period(orbits[name].times, r)
	println("period\t ", T)
	nperi = -orbits[name].times[1] * lguys.T0 / T
	println("nperi\t ", nperi)

	println("f ", orbits[name].f[end])
	println()
end

# ╔═╡ f4de85f0-1583-4d7a-96ed-c4920109600b
begin 
	plots = []
	for name in names
		p[] = plot(orbits[name].times * lguys.T0, orbits[name].rs, label=name)
		push!(plots, p[])
	end
	plots
end

# ╔═╡ 24775da4-4f9e-4545-b414-71c579680d58
R = [
	−0.7500 −0.4572 0.4780 
	−0.2664 0.8702 0.4144 
	0.6054 −0.1835 0.7745
]

# ╔═╡ aeede747-b438-46ec-a11a-15c324112587
function to_stream(ra::Real, dec::Real)
	x = [cosd(ra) * cosd(dec), sind(ra) * cosd(dec), sind(dec)]

	x_t = R * x
	ϕ1 = atand(x_t[2], x_t[1])
	r = sqrt(x_t[1]^2 + x_t[2]^2) # equals cos phi2
	ϕ2 = atand(x_t[3], r)
	return ϕ1, ϕ2
end

# ╔═╡ 83b395c6-bc2d-4146-8bac-d3b708de609b
function to_stream(ra, dec)
	phis = to_stream.(ra, dec)
	phi1 = first.(phis)
	phi2 = last.(phis)
	return phi1, phi2
end

# ╔═╡ da1f56dc-367b-43b0-9d23-7022edb575aa
phi1, phi2 = to_stream(jax.RA, jax.DEC)

# ╔═╡ 1d7d246d-8aab-4c1d-b374-e5d6a79913e2
begin 
	function plot_stream_phi!(obs; kwargs...)
		ra = [o.ra for o in obs]
		dec = [o.dec for o in obs]
		phi1, phi2 = to_stream(ra, dec)

		scatter!(phi1, phi2; kwargs...)
	end
end

# ╔═╡ b236d15b-f07d-44bd-864d-9db325d1d1dc
gr()

# ╔═╡ 03923ddc-c78a-4082-9a79-1f0cf348adb5
pythonplot()

# ╔═╡ 5f4c271d-8912-4c98-8a4c-353a66a7c158
begin 
	p[] = plot(size=(700, 200), xlabel=L"\phi_1\,/\,^{\circ}", ylabel=L"\phi_2\,/\,^{\circ}", legend_position=:bottomright)
	plot_stream_phi!(obss["fiducial"], ms=1, alpha=0.1, label="simulation")
	scatter!(phi1, phi2, ms=1.3, label="Jensen+2021")

	savefig("stream_phi.pdf")
	p[]
end

# ╔═╡ f6592155-99a3-487a-b344-27fb8aaf8e6d
pythonplot()

# ╔═╡ ded33915-5ce9-42de-a492-bcbef96a0700
begin 
	p[] = plot(size=(700, 200), xlabel=L"\phi_1\,/\,^{\circ}", ylabel=L"\phi_2\,/\,^{\circ}", legend_position=:bottomleft)
	plot_stream_phi!(obss["fiducial"], ms=1.5, alpha=0.3, label="simulation", xlims=(-20, 20), ylims=(-5, 5))
	scatter!(phi1, phi2, label="Jensen+2021", shape=:star4)

	savefig("stream_phi_centre.pdf")
	p[]
end

# ╔═╡ 85637c3e-f9d1-446f-b352-1a9b1875cbb0
begin 
	p[] = plot(size=(300, 300), xlabel=L"\phi_1\,/\,^{\circ}", ylabel=L"\phi_2\,/\,^{\circ}", legend_position=:bottomleft)
	plot_stream_phi!(obss["fiducial"], ms=1.5, alpha=0.3, label="simulation", xlims=(-1, 1), ylims=(-1, 1), aspect_ratio=1)

	savefig("very_centre.pdf")
	p[]
end

# ╔═╡ a5849a9e-8105-4c9a-a332-e2ad88fa4dfb
begin 
	plot(size=(700, 200), legend=false)
	for name in names
		plot_stream_phi!(obss[name], ms=1, alpha=0.1)
	end
	scatter!(phi1, phi2)
end

# ╔═╡ eaabe215-ea05-4946-aa3a-0266776be143
phi1_p, phi2_p = to_stream()

# ╔═╡ 8895569b-adbd-4d1f-8032-39204054ffcb
scatter(phi1, phi2, aspect_ratio=1, size=(700, 300))

# ╔═╡ 9cab3d7a-4ed8-4b3a-ac04-cd0d8202dc7b
obs = obss["fiducial"]

# ╔═╡ 536afe66-c0e1-491f-9c12-87505c1ed279
begin 
	scatter([o.pm_ra for o in obs], [o.pm_dec for o in obs], ms=1, alpha=0.1)
	scatter!(jax.pmra, jax.pmdec)
end

# ╔═╡ Cell order:
# ╠═db8a38d4-fae2-11ee-18e3-a189bdef1759
# ╠═04b79ff4-1c47-44bc-a7d1-e3e119dccd99
# ╠═681dc280-e68a-4e52-990b-88089e1803e4
# ╠═1bdea423-3869-4746-9655-ac7c43aaf652
# ╠═d6166e1d-501b-4744-b392-e389b02d7640
# ╠═a499624d-3e0e-445d-a4c8-69be93205d42
# ╠═ac0b3eac-015d-4a96-919e-fc90ddb575bf
# ╠═f9b3062d-cae6-4eaf-b9cd-6dc5a79f14ce
# ╠═15ff24b6-533d-454d-acf8-b3e0f25fed39
# ╠═e7429c5e-216b-4bf6-8dd7-61754af54523
# ╠═8c62cdef-4b48-4d96-add6-05fc6e370372
# ╠═8ff01570-4508-49db-869d-68fab24b3fc0
# ╠═8104512d-c99c-433a-8afc-db949740f3ad
# ╠═5136a85e-5f19-4604-b1d1-5ef7e04bdb29
# ╠═6ae8ec6c-dc20-497e-9fd2-a8e4cffd920a
# ╠═97d8f567-bcbe-4d24-a609-3f2002b71b29
# ╠═0d2c4538-16e7-41d3-b8d8-bd4872f15026
# ╠═f994c87e-15df-421d-9282-361d1bb55f16
# ╠═124e0435-ccad-4918-bebc-762fa846ee69
# ╠═de6d5e1e-2937-4199-9274-a0bc49fb8d8a
# ╠═dbc15de6-8c9b-4e25-92dc-c02a68976355
# ╠═4e103842-f084-4a53-a5ae-c149447744a1
# ╠═1fadb10f-b160-41f7-99b3-25e654d0a674
# ╠═b90632da-680c-49d5-9903-b00191a4f03b
# ╠═ed11b804-49f8-463b-84fd-7b3c18a0c5f6
# ╠═3813d261-6890-44aa-9560-c2758745eec4
# ╠═b39e1baf-2de0-4539-9957-200679fcd38d
# ╠═a3acf540-ac34-49f3-8594-6d8c4718c83e
# ╠═38d1adb5-a767-429f-bf6a-644e23bb0655
# ╠═fcbbfb05-f761-42c2-a6e6-1121687b0a62
# ╠═874346e2-d285-42b7-aa94-c34d139cb2ee
# ╠═db5ad747-2e70-4861-a260-7bf0496eefd9
# ╠═c30590cb-1d28-438d-8c98-31e752303fee
# ╠═4d71c42f-7bfe-436b-aae0-9ffe75cbfc81
# ╠═8ae6acb2-c4fa-4718-9b44-1937d04be1f0
# ╠═368483d0-84d4-48ef-8f61-95ac001a587f
# ╠═6807c973-c193-4f66-9a1e-a569f3ad872f
# ╠═6fecdb3a-9945-4027-ba21-07bd6fec9a36
# ╠═c649a8a6-2ede-4988-ba69-c6df360931a5
# ╠═cf2f33a0-d940-4ddb-bc8c-592af2dff0a7
# ╠═6de087be-fa5a-4f2e-ab5d-189385dc2676
# ╠═41b5f8c4-f514-4c83-a46e-2c4e7fc8b492
# ╠═cc81994f-3944-4e01-bb8e-cfc2982c800c
# ╠═3d02ec0a-7376-4530-a66f-7297454d98e1
# ╠═dc7c775a-606b-4fb4-a688-4c858c31ca20
# ╠═d36a3b1c-d224-4960-9fa6-aab119b8ab61
# ╠═f201f6c4-f2dc-4e28-91b7-cd3f5ce871b2
# ╠═05bd3bd4-93dc-48d8-ac97-ecddc9390fca
# ╠═f4de85f0-1583-4d7a-96ed-c4920109600b
# ╠═24775da4-4f9e-4545-b414-71c579680d58
# ╠═aeede747-b438-46ec-a11a-15c324112587
# ╠═83b395c6-bc2d-4146-8bac-d3b708de609b
# ╠═da1f56dc-367b-43b0-9d23-7022edb575aa
# ╠═1d7d246d-8aab-4c1d-b374-e5d6a79913e2
# ╠═b236d15b-f07d-44bd-864d-9db325d1d1dc
# ╠═03923ddc-c78a-4082-9a79-1f0cf348adb5
# ╠═5f4c271d-8912-4c98-8a4c-353a66a7c158
# ╠═f6592155-99a3-487a-b344-27fb8aaf8e6d
# ╠═ded33915-5ce9-42de-a492-bcbef96a0700
# ╠═85637c3e-f9d1-446f-b352-1a9b1875cbb0
# ╠═a5849a9e-8105-4c9a-a332-e2ad88fa4dfb
# ╠═eaabe215-ea05-4946-aa3a-0266776be143
# ╠═8895569b-adbd-4d1f-8032-39204054ffcb
# ╠═9cab3d7a-4ed8-4b3a-ac04-cd0d8202dc7b
# ╠═536afe66-c0e1-491f-9c12-87505c1ed279
