### A Pluto.jl notebook ###
# v0.20.5

using Markdown
using InteractiveUtils

# ╔═╡ 56f02c6e-0a60-11f0-37e4-db5ba2ec1ba4
begin 
	using Pkg; Pkg.activate()

	using LilGuys
end

# ╔═╡ 2c61ebbe-b228-4330-9f80-313afd61f28f
using CairoMakie, Arya

# ╔═╡ 8f5a0186-c8ee-48fc-89e7-1637f5c44280
function Ms_from_Mv(Mv, Y)
	Ls = 10 ^ (-0.4 * (Mv - 4.83)) * Y
end

# ╔═╡ 8151de8a-3c1e-4010-beae-e33607f43fba
Ms_from_Mv(-10.82, 1.5)

# ╔═╡ 24ef95d2-78c1-40bf-a9e0-df34df4e2ae9
function halo_from_Mv(Mv, Y=1.5, δr=0)
	v_max = LilGuys.vel_from_M_s_fattahi(Ms_from_Mv(Mv, Y)/M2MSUN)
	r_max = LilGuys.Ludlow.solve_rmax(v_max)
	return LilGuys.TruncNFW(r_circ_max=r_max*10^δr, v_circ_max=v_max, trunc=10)
end

# ╔═╡ b392dd40-666f-4d50-aae0-834865961bfb
halo = halo_from_Mv(-10.82, 1.5)

# ╔═╡ 72fc5d32-8f00-4a58-a81a-43099d5affcd
LilGuys.mass(halo, 0.1) / LilGuys.mass(halo, 20)

# ╔═╡ 07c8b513-7344-4267-b25e-993137a6996e
α = LilGuys.R_h(LilGuys.Exp2D())

# ╔═╡ 9caf94ba-34eb-4e23-8346-ed5d230c894a
function M_h_obs(halo, R_h)
	
	σ = LilGuys.σv_star_mean(halo, LilGuys.Exp2D(R_s=R_h/α)) 

	return 5/2 * σ^2 * R_h  / LilGuys.G
end

# ╔═╡ 9d4feebe-c728-451e-b9af-8b8d6cc9742c
LilGuys.σv_star_mean(halo, LilGuys.Exp2D(R_s=0.1/α)) 

# ╔═╡ c33b100e-ce25-423e-8bd2-aa972f055558
sqrt(LilGuys.v_circ(halo, 0.3) * 0.3)

# ╔═╡ da79ba87-a73d-4cd5-9cdf-99ae23921c2a
M_h_obs(halo,0.3)

# ╔═╡ 814619b7-9d4c-48cf-9e5a-fab6fcd972a7
mass(halo, 0.3)

# ╔═╡ 2c3fa7fd-1334-443a-8794-bfafcd648e33
CairoMakie.activate!(type=:png)

# ╔═╡ 40fc5d6b-a8a7-472f-bdab-b285d287bfe1
log10(2)

# ╔═╡ 7b1692af-74b0-447d-b138-d73df8bb6e0a
function log_mass_ratio(log_R_h, Mv, δr=0)
	R_h = 10 .^ log_R_h

	halos = halo_from_Mv(Mv, 1.5, δr)
	ms_o = M_h_obs(halos, R_h)
	ms_r = mass(halos, 30)
	return log10(ms_o/ms_r)
end

# ╔═╡ 205daee5-c379-47db-b977-641deabb8c3d
function log_mass_ratio_h(log_R_h, Mv, δr=0)
	R_h = 10 .^ log_R_h

	halos = halo_from_Mv(Mv, 1.5, δr)
	ms_o = M_h_obs(halos, R_h)
	ms_r = mass(halos, R_h)
	return log10(ms_o/ms_r)
end

# ╔═╡ 12730072-41f0-4a9b-9831-f69ea15ff790
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = L"\log\,R_h",
		ylabel = L"\log\,M_\textrm{est} / M_\textrm{tot}"
	)

	N = 20
	log_R_h = LinRange(-3, 2, N)
	Mvs = [-10, -8, -6, -4, -2, 2]

	colorrange = extrema(Mvs)
	for Mv in Mvs
		lines!(log_R_h, x->log_mass_ratio(x,Mv) , 
				   linestyle=:solid, color=Mv, colorrange=colorrange,
				   label = "$Mv"
				  )

	end

	axislegend("Mv", position=:rb)
	
	fig
end

# ╔═╡ 550f081d-12c6-45d8-8b12-3a2714f899b6
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = L"\log\, R_h",
		ylabel = L"\log M_\textrm{est} / M_\textrm{h}"
	)

	N = 20
	log_R_h = LinRange(-3, 2, N)
	Mvs = [-10, -8, -6, -4, -2, 2]

	colorrange = extrema(Mvs)
	for Mv in Mvs
		lines!(log_R_h, x->log_mass_ratio_h(x,Mv) , 
				   linestyle=:solid, color=Mv, colorrange=colorrange,
				   label = L"%$Mv"
				  )

	end

	axislegend("Mv", position=:lb)
	
	fig
end

# ╔═╡ 9810993a-9e8e-4a34-a09a-c6eb56bdeef1
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log R h",
		ylabel = L"\log M_\textrm{est} / M_\textrm{tot}"
	)

	N = 20
	log_R_h = LinRange(-2, 1, N)

	Mv = -8

	δrs = LinRange(-0.2, 0.2, 5)
	colorrange = extrema(δrs)

	for δr in δrs
		lines!(log_R_h, x->log_mass_ratio(x,Mv, δr), linestyle=:solid, color=δr, colorrange=colorrange, label=string(round(δr, digits=2)))
	end

	axislegend("δr", position=:lt)

	fig
end

# ╔═╡ 2ae77743-8a8e-4913-9f55-4ee7d338d26a
let
	fig = Figure()
	ax = Axis(fig[1,1],
		xlabel = "log R h",
		ylabel = L"\log M_\textrm{est} / M_\textrm{h}"
	)

	N = 20
	log_R_h = LinRange(-2, 1, N)

	Mv = -8

	δrs = LinRange(-0.2, 0.2, 10)
	colorrange = extrema(δrs)

	for δr in δrs
		lines!(log_R_h, x->log_mass_ratio_h(x,Mv, δr), linestyle=:solid, color=δr, colorrange=colorrange)
	end

	fig
end

# ╔═╡ Cell order:
# ╠═56f02c6e-0a60-11f0-37e4-db5ba2ec1ba4
# ╠═2c61ebbe-b228-4330-9f80-313afd61f28f
# ╠═8f5a0186-c8ee-48fc-89e7-1637f5c44280
# ╠═8151de8a-3c1e-4010-beae-e33607f43fba
# ╠═b392dd40-666f-4d50-aae0-834865961bfb
# ╠═24ef95d2-78c1-40bf-a9e0-df34df4e2ae9
# ╠═72fc5d32-8f00-4a58-a81a-43099d5affcd
# ╠═07c8b513-7344-4267-b25e-993137a6996e
# ╠═9caf94ba-34eb-4e23-8346-ed5d230c894a
# ╠═9d4feebe-c728-451e-b9af-8b8d6cc9742c
# ╠═c33b100e-ce25-423e-8bd2-aa972f055558
# ╠═da79ba87-a73d-4cd5-9cdf-99ae23921c2a
# ╠═814619b7-9d4c-48cf-9e5a-fab6fcd972a7
# ╠═2c3fa7fd-1334-443a-8794-bfafcd648e33
# ╠═40fc5d6b-a8a7-472f-bdab-b285d287bfe1
# ╠═7b1692af-74b0-447d-b138-d73df8bb6e0a
# ╠═205daee5-c379-47db-b977-641deabb8c3d
# ╠═12730072-41f0-4a9b-9831-f69ea15ff790
# ╠═550f081d-12c6-45d8-8b12-3a2714f899b6
# ╠═9810993a-9e8e-4a34-a09a-c6eb56bdeef1
# ╠═2ae77743-8a8e-4913-9f55-4ee7d338d26a
