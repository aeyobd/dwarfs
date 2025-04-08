### A Pluto.jl notebook ###
# v0.20.5

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

# ╔═╡ 25f6474b-2ee7-4d2e-9ce5-fdf782adfa72
using Agama

# ╔═╡ dea35f32-ba19-474f-a7b9-fdaa5e89fc43
include("./paper_style.jl")

# ╔═╡ 2a7db353-9975-41d0-8a03-d91622d549d4
potential_dir = joinpath(ENV["DWARFS_ROOT"], "agama", "potentials")

# ╔═╡ 3ef82ecd-0eba-4905-b370-411ea2f2384b
readdir(potential_dir)

# ╔═╡ 5cbabcc8-6473-4ac1-8084-80da24a84617
Φ_ep = Agama.Potential(file = joinpath(potential_dir, "EP2020.ini"))

# ╔═╡ 44fa13a8-0119-49f9-8430-eb46d07430ca
Φ_mm = Agama.Potential(file = joinpath(potential_dir, "mcmillan17.ini"))

# ╔═╡ 8cc0f250-f7b0-4336-9d4d-d9fa24217b74
Φ_lmc = Agama.Potential(file = joinpath(potential_dir, "vasiliev24/L3M11/potential_lmc_init.ini"))

# ╔═╡ aca16e3c-b1bc-42ee-84b4-0fba6c5ebbf3
Φ_v24 = Agama.Potential(file = joinpath(potential_dir, "vasiliev24/L3M11/potential_mw_now.ini"))

# ╔═╡ babbbe3b-13b7-4450-98be-f807f21706d8
Φ_bulge, Φ_thin, Φ_thick, Φ_nfw = Φ_ep._py

# ╔═╡ dc9ac37f-0fdf-4061-b924-af77ecd21396
module AU 
	include(joinpath(ENV["DWARFS_ROOT"], "utils", "agama_utils.jl"))

end

# ╔═╡ e01b33e6-8c15-421c-a1ad-d5f6fbc3997a
AU.V_V2KMS

# ╔═╡ 86961f65-fa3e-4899-81d5-94c6bec3b8de
function plot_v_circ!(pot; x=LinRange(-1, 3, 1000), kwargs...)
	r = 10 .^ x

	m = pot.enclosedMass(r) |> Agama.py2vec
	v = @. sqrt( m / r) * V2KMS

	lines!(x, v; kwargs...)
	
end

# ╔═╡ 9a81c578-b5e9-45de-9ee5-5c81d52d8871
function plot_v_circ_v!(pot; x=LinRange(-1, 3, 1000), kwargs...)
	r = 10 .^ x ./ AU.V_R2KPC

	m = pot.enclosedMass(r) |> Agama.py2vec
	v = @. sqrt( m / r) * AU.V_V2KMS

	lines!(x, v; kwargs...)
	
end

# ╔═╡ cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
CairoMakie.activate!(type=:png)

# ╔═╡ 3c032178-8d48-4f9c-bcec-9bf704718ea9
@savefig "v_circ_potential" let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$\log\, R$ / kpc",
		ylabel = L"$v_\textrm{circ}$ / km\,s$^{-1}$",
		yscale=log10, 
		yticks = [10, 20, 100, 200],
		yminorticks = [10:10:100; 100:100:1000;]
	)

	plot_v_circ!(Φ_ep._py, label="total", linestyle=:solid, color=:black, linewidth=1.5)

	plot_v_circ!(Φ_bulge, label="bulge", linestyle=:dot, color=COLORS[1])
	plot_v_circ!(Φ_thin + Φ_thick, label = "disk", linestyle=:dashdot,  color=COLORS[2])
	plot_v_circ!(Φ_nfw, label="halo", linestyle=:dash,  color=COLORS[3])


	plot_v_circ_v!(Φ_v24._py, label="V24", linestyle=:solid, color=COLORS[4])
	plot_v_circ_v!(Φ_lmc._py, label="LMC", color=COLORS[5])

	axislegend(position=:cb)
	#Legend(fig[1,2], ax)
	fig
end

# ╔═╡ 471d64b6-0176-4db6-8603-9de191997b9d
let
	fig = Figure()

	ax = Axis(fig[1,1], 
		xlabel = L"$\log\, R$ / kpc",
		ylabel = L"$v_\textrm{circ}$ / km\,s$^{-1}$",
		yscale=log10, 
		yticks = [30, 100, 300],
		yminorticks = [10:10:100; 100:100:1000;]
	)

	plot_v_circ!(Φ_ep._py, label="total", linestyle=:solid, color=:black, linewidth=1.5)

	plot_v_circ_v!(Φ_v24._py, label="V24")
	plot_v_circ_v!(Φ_lmc._py, label="LMC")
	plot_v_circ!(Φ_mm._py, label="MM")


	axislegend(position=:lb)
	fig
end

# ╔═╡ Cell order:
# ╠═0125bdd2-f9db-11ef-3d22-63d25909a69a
# ╠═25f6474b-2ee7-4d2e-9ce5-fdf782adfa72
# ╠═dea35f32-ba19-474f-a7b9-fdaa5e89fc43
# ╠═2a7db353-9975-41d0-8a03-d91622d549d4
# ╠═3ef82ecd-0eba-4905-b370-411ea2f2384b
# ╠═5cbabcc8-6473-4ac1-8084-80da24a84617
# ╠═44fa13a8-0119-49f9-8430-eb46d07430ca
# ╠═8cc0f250-f7b0-4336-9d4d-d9fa24217b74
# ╠═aca16e3c-b1bc-42ee-84b4-0fba6c5ebbf3
# ╠═babbbe3b-13b7-4450-98be-f807f21706d8
# ╠═dc9ac37f-0fdf-4061-b924-af77ecd21396
# ╠═e01b33e6-8c15-421c-a1ad-d5f6fbc3997a
# ╠═86961f65-fa3e-4899-81d5-94c6bec3b8de
# ╠═9a81c578-b5e9-45de-9ee5-5c81d52d8871
# ╠═cf74dd35-70f3-4a88-b1c0-580d0ffe69aa
# ╠═3c032178-8d48-4f9c-bcec-9bf704718ea9
# ╠═471d64b6-0176-4db6-8603-9de191997b9d
