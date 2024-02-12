### A Pluto.jl notebook ###
# v0.19.36

using Markdown
using InteractiveUtils

# ╔═╡ bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
begin 
	using Pkg; Pkg.activate()
	import LilGuys as lguys
	
	using Plots; plotly()
end

# ╔═╡ 6964cbd2-9ef9-4682-9ac1-80c10518c374
out = lguys.Output("out")

# ╔═╡ a558298a-a698-4522-9f48-1a6a0b074814
snap = out[1000]

# ╔═╡ 3f0bad35-88dc-4fb7-90c4-1b5395ab53a4
function scatter_snap(points...; kwargs...)
	plots = lguys.scatter_xyz(snap.positions, msw=0, ms=1, ma=0.1, mc="black")

	for i in 1:3
		for point in points
			
			scatter([point[i]], [point[(i+1) % 4 + 1]]; kwargs...)
		end
	end
	return plots
end

# ╔═╡ 6e42d56b-82e2-4974-bdfc-cdd07c0874f6
md"""
# The density centre

The most direct algorithm is to find the point of highest density. However, we do have to estimate a continuous physical quantity from 
"""

# ╔═╡ 0ffd8ac7-6ee3-47f4-a83a-59f9bdb3ee85
c_p = lguys.centre_ρ(snap)

# ╔═╡ 0de3cc69-6c25-4a8b-bfa4-34987e0da961
 ρs = lguys.calc_ρs(snap, η=1/√5, k=5)

# ╔═╡ 085fc1e3-ed40-45eb-91fc-7ca4275637c7
c_i_ρ =  snap.positions[:, argmax(ρs)]

# ╔═╡ 281b2b17-55f6-4750-94b0-eff1cbe7ca56
histogram(log10.(ρs))

# ╔═╡ 7ef52925-8ada-4e2d-ad2a-679de9d23ca9
lguys.scatter_xyz(snap.positions, msw=0, ms=1, ma=0.1, mc="black")

# ╔═╡ 89f3e4cc-3bc5-4b6e-a88d-be8f84e35ad5
begin scatter_
	scatter(snap.positions[2, :], snap.positions[3, :], marker_z=log10.(ρs), msw=0, ms=1)
	scatter!([c_p.minimizer[2]], [c_p.minimizer[3]])
	scatter!([c_i_ρ[2]], [c_i_ρ[3]])

end

# ╔═╡ 6e9a64c9-dd0d-4044-87a6-ee412b5d6f39
scatter_snap(c_p.minimizer)

# ╔═╡ 6d623271-3de6-4527-9b43-1eb5546489d2
begin 
	scatter(snap.positions[1, :], snap.positions[3, :], marker_z=log10.(ρs), msw=0, ms=1)
	scatter!([c_p.minimizer[1]], [c_p.minimizer[3]])
	scatter!([c_i_ρ[1]], [c_i_ρ[3]])

end

# ╔═╡ 13a9e76a-f6f6-4bf0-aa43-a485c95e9fa2
md"""
# The potential centre
Another simple one, which point has the minimum potential? 

"""

# ╔═╡ 4431ea6c-84b2-416d-90b2-5cc238a0c841
md"""
# The shrinking spheres centre

Algorithm

Repeat
- Find centroid
- Remove particles at some percentile of radii
- Remove unbound particles according to the COM
"""

# ╔═╡ b37cf8a0-1125-42ac-9482-57516b4d3b19
c_p.minimizer

# ╔═╡ 68a6f529-bb98-41f4-b2d6-b14dc4b4a2b0
c_i_ρ

# ╔═╡ 693bcb86-29b8-4dd6-989b-641af3ff4b42
p = Ref{Plots.Plot}()

# ╔═╡ 6fcd5038-3337-4f48-bb36-5f7cc07334a0
function ϕ_eff(positions, velocities, masses; h=0.08)
	N = length(masses)
	ϕ = zeros(N)

	for i in 1:N
		if i % 100 == 0
			print("\r $i / $N")
		end

		ϕ_g = -masses ./ sqrt.(lguys.calc_r(positions[:, i], positions).^2 .+ h.^2)
		ϕ_v =  1/N * 1/2 * lguys.calc_r(velocities[:, i], velocities).^2
			
		ϕ[i] = sum(min.(ϕ_g .+ ϕ_v, 0))
		
	end
	return ϕ
end

# ╔═╡ 0ef2998e-1207-4d44-a725-367f4eab010f
ps = ϕ_eff(snap.positions, snap.velocities, snap.masses)

# ╔═╡ 20624c12-b10b-43d7-9afe-df9518522172
ps2 = ϕ_eff(snap.positions, snap.velocities, -ps)

# ╔═╡ 2be8f04f-79ec-4824-a148-636a70913e8a
cutoff = lguys.percentile(ps, 10)

# ╔═╡ bd8444b2-5c7e-47f4-a3c8-e504038f82ae


# ╔═╡ 51044823-4aa9-4a54-b3bb-cc6ef6d98ee5
ss_c = lguys.ss_centre(snap, history=true)

# ╔═╡ df5d1b78-2648-473c-9ce3-7c593d685313
phi_filt = ps .< cutoff

# ╔═╡ b6026809-a596-437e-bf4a-e74b6aa88a32
ss_xs = hcat(ss_c.history.xs...)

# ╔═╡ 7ec18905-02a7-4142-9346-8f1e8882ae16
x_c = lguys.centroid(snap.positions[:, phi_filt], (ps[phi_filt] .- cutoff).^2)

# ╔═╡ da1330c4-5e04-4cb4-85d6-594012672a70
begin 
	plot()
	scatter!(snap.positions[2, :], snap.positions[3, :], ms=1, color="black", alpha=0.05)
		scatter!(ss_xs[2, :], ss_xs[3, :], marker_z=1:size(ss_xs, 2), msw=0)

	scatter!(x_c[2, :], x_c[3, :])
	xlims!(-200, 200)
	ylims!(-200, 200)
	# savefig("centre_finding")
end

# ╔═╡ 6d7f69be-1f9e-4f8a-8131-6618875f5581
begin 
	scatter(snap.positions[3, phi_filt], snap.positions[2, phi_filt], marker_z=ps[phi_filt], ms=3, msw=0, ma=1)
	scatter!(x_c[3, :], x_c[2, :])
end

# ╔═╡ bc259579-8cfe-4646-8856-6febec46f916
scatter(snap.positions[3, phi_filt], snap.positions[2, phi_filt], marker_z=snap.Φs[phi_filt], ms=2, msw=0, ma=10.1)

# ╔═╡ 05fb9044-9704-4f47-98ef-61b98a013ffd
begin 
	scatter(snap.positions[3, :], snap.positions[2, :], marker_z=snap.Φs, ms=2, msw=0, ma=1)
	xlims!(-200, 200)
	ylims!(-200, 200)
end

# ╔═╡ 8089c884-753c-4ba6-8f21-2793e79a2ef4
begin 
	scatter(snap.positions[2, :], snap.positions[3, :], marker_z=ps, ms=2, msw=0, ma=1)
		xlims!(-200, 200)
	ylims!(-200, 200)
end

# ╔═╡ b01611ca-040e-42be-8a65-4b8eb6b8a63c
begin 
	histogram(ps, yscale=:log10)
	vline!([cutoff, 0])
end

# ╔═╡ 97dc0483-8841-4fc8-be0b-1b84c54d74b7
histogram(lguys.calc_r(snap.velocities))

# ╔═╡ 01a3e72b-fcdc-4d28-b5ce-2598ea457f94
histogram(snap.Φs)

# ╔═╡ 9a69f8ab-2f5b-4cfa-9056-e36bb5bfc954
minimum(ps)

# ╔═╡ Cell order:
# ╠═bb92b6c2-bf8d-11ee-13fb-770bf04d91e9
# ╠═6964cbd2-9ef9-4682-9ac1-80c10518c374
# ╠═a558298a-a698-4522-9f48-1a6a0b074814
# ╠═3f0bad35-88dc-4fb7-90c4-1b5395ab53a4
# ╟─6e42d56b-82e2-4974-bdfc-cdd07c0874f6
# ╠═0ffd8ac7-6ee3-47f4-a83a-59f9bdb3ee85
# ╠═085fc1e3-ed40-45eb-91fc-7ca4275637c7
# ╠═0de3cc69-6c25-4a8b-bfa4-34987e0da961
# ╠═281b2b17-55f6-4750-94b0-eff1cbe7ca56
# ╠═7ef52925-8ada-4e2d-ad2a-679de9d23ca9
# ╠═89f3e4cc-3bc5-4b6e-a88d-be8f84e35ad5
# ╠═6e9a64c9-dd0d-4044-87a6-ee412b5d6f39
# ╠═6d623271-3de6-4527-9b43-1eb5546489d2
# ╟─13a9e76a-f6f6-4bf0-aa43-a485c95e9fa2
# ╟─4431ea6c-84b2-416d-90b2-5cc238a0c841
# ╠═b37cf8a0-1125-42ac-9482-57516b4d3b19
# ╠═68a6f529-bb98-41f4-b2d6-b14dc4b4a2b0
# ╠═693bcb86-29b8-4dd6-989b-641af3ff4b42
# ╠═6fcd5038-3337-4f48-bb36-5f7cc07334a0
# ╠═0ef2998e-1207-4d44-a725-367f4eab010f
# ╟─20624c12-b10b-43d7-9afe-df9518522172
# ╠═2be8f04f-79ec-4824-a148-636a70913e8a
# ╠═bd8444b2-5c7e-47f4-a3c8-e504038f82ae
# ╠═51044823-4aa9-4a54-b3bb-cc6ef6d98ee5
# ╠═df5d1b78-2648-473c-9ce3-7c593d685313
# ╠═b6026809-a596-437e-bf4a-e74b6aa88a32
# ╠═da1330c4-5e04-4cb4-85d6-594012672a70
# ╠═7ec18905-02a7-4142-9346-8f1e8882ae16
# ╠═6d7f69be-1f9e-4f8a-8131-6618875f5581
# ╠═bc259579-8cfe-4646-8856-6febec46f916
# ╠═05fb9044-9704-4f47-98ef-61b98a013ffd
# ╠═8089c884-753c-4ba6-8f21-2793e79a2ef4
# ╠═b01611ca-040e-42be-8a65-4b8eb6b8a63c
# ╠═97dc0483-8841-4fc8-be0b-1b84c54d74b7
# ╠═01a3e72b-fcdc-4d28-b5ce-2598ea457f94
# ╠═9a69f8ab-2f5b-4cfa-9056-e36bb5bfc954
