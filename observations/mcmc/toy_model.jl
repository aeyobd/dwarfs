### A Pluto.jl notebook ###
# v0.20.18

using Markdown
using InteractiveUtils

# ╔═╡ 1368181a-97d7-11f0-3b4f-b90907cb3133
begin
	import Pkg; Pkg.activate()

	using LilGuys
	using Arya
	using DataFrames, CSV
	import TOML
	using CairoMakie
end

# ╔═╡ 19112874-f462-4cfb-bbef-ca2a43e72fd1
using StatsBase, Distributions

# ╔═╡ 7288a3ad-375a-4930-a98b-2d225d8ed709
using PyFITS

# ╔═╡ f1b3eef7-7aae-4622-ba01-9470a2867b39
f_mix = 1

# ╔═╡ 446bb406-8ec0-4dba-822b-187c0b97d73b
N_stars = 100_000

# ╔═╡ cc93ac36-795d-443e-998e-968a2d01f858
R_max = 30

# ╔═╡ a235a388-a42c-494f-8ef9-26554c5238b8
f_sat_true = 0.03

# ╔═╡ 82ff3270-fbe2-45df-943c-946f1292075a
P_cut = 0.5

# ╔═╡ 688e1152-81db-4a93-a314-1d902c87f054
x_dist_memb = Beta(0.1, 8)

# ╔═╡ 0feed8ce-18c4-4213-8519-5749a4b16552
prof_memb_assumed = LilGuys.Plummer(r_s=3.0)

# ╔═╡ ce09dc9f-0cea-4121-b0da-60a808399b2f
prof_memb = LilGuys.Exp2D(R_s=0.8)

# ╔═╡ 86671511-ac67-4837-bc29-1f4c8597519f
N_memb = ceil(Int, N_stars*f_sat_true)

# ╔═╡ a6af0647-c18e-4157-bf71-ed63186fe546
N_bkd = ceil(Int, N_stars*(1-f_sat_true))

# ╔═╡ 2ee203cd-1370-441f-96de-bf086ca48bc6
f_sat = f_sat_true

# ╔═╡ b719c6db-ec7c-45c3-a1d1-da9778c47477
md"""
# Toy Model
"""

# ╔═╡ 7d78cfc8-18b9-416b-9676-1c86c9f84841
df_scl = let 
	df = read_fits("../jax_data/sculptor.fits")

	df[!, :PSAT_NO] = df.L_PM_SAT .* df.L_CMD_SAT ./ (df.L_PM_SAT .* df.L_CMD_SAT .+ df.L_PM_BKD .* df.L_CMD_BKD)

	df[df.F_BEST .=== 1.0, :]
end

# ╔═╡ b650ec10-1101-4032-8f42-34fb2b9681d2
hist(log10.(df_scl.L_PM_SAT .* df_scl.L_CMD_SAT ./df_scl.L_PM_BKD .* df_scl.L_CMD_BKD)[df_scl.PSAT .> 0], bins=LinRange(-5, 5, 100))

# ╔═╡ 6ae2719c-f29a-4e6a-875e-90eed9c3a998
logit(x) = log(x / (1-x))

# ╔═╡ f2b73277-eccb-4e6c-8ec4-a76533e9fcdb
L_dist_memb = x -> pdf(x_dist_memb, x)

# ╔═╡ 8e864e73-1c7f-4139-b921-db18fa004120
LilGuys.integrate(x -> (L_dist_memb(x) / (1 + L_dist_memb(x)) > P_cut), 0, 1) * f_sat

# ╔═╡ 99b70a40-b66b-4f24-b66e-f8e09aeca0f1
L_dist_memb

# ╔═╡ 7e827a7e-99fe-452d-be3e-24ba1a741887
df_bkd = let
	df = DataFrame(
		:R => sqrt.(rand(Uniform(0, R_max ^ 2), N_bkd))
	)
	θ_bkd = rand(Uniform(0, 2π), N_bkd)

	df[!, :x] = df.R .* cos.(θ_bkd)
	df[!, :y] = df.R .* sin.(θ_bkd)

	df[!, :x_nospace] = rand(Uniform(0, 1), N_bkd)

	df[!, :L_memb_nospace] = L_dist_memb(df.x_nospace)
	df
end

# ╔═╡ f73aa823-08a9-4e21-bede-a9095647cd5a
Σ_memb = R -> LilGuys.surface_density(prof_memb, R)

# ╔═╡ 0afcd462-ab26-4b0c-a27a-3a838b0225e3
df_memb = let
	df = DataFrame(
		:R => LilGuys.sample_surface_density(Σ_memb, N_memb),
		:x_nospace => (rand(x_dist_memb, N_memb)),
			#rand(N_memb) .< completeness
	)

	df[!, :L_memb_nospace] = L_dist_memb(df.x_nospace)
	θ = rand(Uniform(0, 2π), N_memb)
	df[!, :x] = df.R .* cos.(θ)
	df[!, :y] = df.R .* sin.(θ)

	df
end


# ╔═╡ 5ed34292-2f95-4730-abcd-cb60a5509982
hist( df_memb.L_memb_nospace ./ (df_memb.L_memb_nospace .+ 1))

# ╔═╡ db089c17-a106-4b61-9cff-8c6c04971355
let
	fig = Figure()
	ax = Axis(fig[1,1])
	bins = LinRange(-5, 5, 100)

	x = df_scl.PSAT_NO[df_scl.PSAT .> P_cut]
	hist!(logit.(x), normalization=:pdf, bins=bins)

	x = df_memb.L_memb_nospace ./ (df_memb.L_memb_nospace .+ 1)
	stephist!(logit.(x), normalization=:pdf, bins=bins, color=COLORS[3])



	ax = Axis(fig[2,1])


	x = df_scl.PSAT_NO[0 .< df_scl.PSAT .< P_cut]
	hist!(logit.(x), normalization=:pdf, bins=bins, color=COLORS[2])

	x = df_bkd.L_memb_nospace ./ (df_bkd.L_memb_nospace .+ 1)
	stephist!(logit.(x), normalization=:pdf, bins=bins, color=COLORS[4])

	
	fig
end

# ╔═╡ d72a5173-5b76-4f30-8578-e615a5efecff
hist(df_memb.L_memb_nospace)

# ╔═╡ abfce12c-a3b2-4699-89db-734c87f1dd85
scatter(df_bkd.x, df_bkd.y, markersize=1)

# ╔═╡ 59466065-3e1b-4cf2-9649-c48a10dfed05
scatter(df_memb.x, df_memb.y, markersize=1)

# ╔═╡ 4fa3e12a-3407-4fcb-8452-b9cc46182d7f
hist(log10.(df_bkd.L_memb_nospace))

# ╔═╡ 8021eb9e-9064-47d8-a108-afafa5c60ac4
hist(log10.(df_memb.L_memb_nospace))

# ╔═╡ 3f327a21-2532-4768-8fbb-4a68f2949a5f
md"""
# Combined sample
"""

# ╔═╡ 6aae7d18-eb3b-4ff8-96a1-06ef7a2ebfd8
Σ_bkd_model = R -> 1 / (π * R_max^2)

# ╔═╡ bf609cbf-4755-4f43-9c8d-9f3601a29cf4
Σ_memb_model = R -> (1-f_mix) *  LilGuys.surface_density(prof_memb, R) .+ f_mix * LilGuys.surface_density(prof_memb_assumed, R)

# ╔═╡ 03a7b639-71a9-4c7f-9d97-57b07abbb3ae
df_all = let
	df = vcat(df_memb, df_bkd)
	df[!, :is_true_memb] = [trues(N_memb); falses(N_bkd)]

	df[!, :L_memb] = Σ_memb_model.(df.R) .* df.L_memb_nospace
	df[!, :L_bkd] =  Σ_bkd_model.(df.R)

	df[!, :P_sat] = @. (f_sat * df.L_memb) / (f_sat * df.L_memb + (1-f_sat) * df.L_bkd)
	df[!, :P_sat_nospace] = @. (f_sat * df.L_memb_nospace) / (f_sat * df.L_memb_nospace + (1-f_sat) * 1)
	df
end

# ╔═╡ 12711fe9-6518-455d-93e2-049d6cda0566
scatter(df_all.x, df_all.y, alpha=0.1, markersize=1)

# ╔═╡ 5fc1c365-2634-4781-8729-7c4c774fb8cc
hist(log10.(df_all.L_bkd ./ df_all.L_memb))

# ╔═╡ bbb6cf08-48e4-4aa6-94a7-853c87d71b50
df_recovered = df_all[df_all.P_sat .> P_cut, :]

# ╔═╡ 365909af-8fa6-49d1-8edb-3c7b646792e3
purity = mean(df_recovered.is_true_memb)

# ╔═╡ 34b0555d-605b-44e8-a467-8fb0962a54c4
completeness = sum(df_recovered.is_true_memb) / N_memb

# ╔═╡ 9b918b65-0bf5-4564-8e09-c797adeb0b05
mean(df_memb.L_memb_nospace), completeness

# ╔═╡ d51390ea-f0cb-43f7-900a-539a901502f2
sum(.!df_recovered.is_true_memb) / N_bkd

# ╔═╡ 3064eeaf-df33-40f3-8f61-2aa3405a7878
scatter(df_all.P_sat_nospace, df_all.P_sat)

# ╔═╡ 4eace3ff-20fd-4a73-8e5e-8b9f9f5a34f5
hist(df_all.P_sat_nospace)

# ╔═╡ 2675f07f-4a2b-46f8-a4a1-0fb994b6ef68
hist(df_all.P_sat)

# ╔═╡ a027b51c-aa7f-4158-80e3-c79f2c866e3d
scatter(df_recovered.x, df_recovered.y)

# ╔═╡ d6f5d1cb-2b0b-427e-b0f9-e2ef874a5eba
md"""
# density profiles
"""

# ╔═╡ 10a4e8f4-8df5-49c3-b0f0-86b751936925
bins=-2:0.1:log10(R_max)

# ╔═╡ 4aed5645-e6a7-4935-8810-20e292149e77
prof = LilGuys.SurfaceDensityProfile(df_recovered.R, bins=bins) |> LilGuys.filter_empty_bins

# ╔═╡ 5c40c51f-690d-4d20-80ce-9446c220015b
prof_true = LilGuys.SurfaceDensityProfile(df_memb.R, bins=bins)

# ╔═╡ ed56f4f0-3946-4434-a613-59ef5efa974a
prof_nospace = LilGuys.SurfaceDensityProfile(df_all.R[df_all.P_sat_nospace .> P_cut], bins=bins)

# ╔═╡ 1b8adcf9-2c19-4d5f-bb3c-5f6e77368a64
prof_all = LilGuys.SurfaceDensityProfile(df_all.R,bins=bins)

# ╔═╡ e6b64d58-b9bb-4523-b38c-48625d7dabcb
f_contam = 1 - purity

# ╔═╡ 61d5f5c3-daf7-4bd8-ac88-69a5ab453ab1
let
	fig = Figure()
	ax = Axis(fig[1,1])
	
	lines!(prof_true.log_R, prof_true.log_Sigma, label="true")
	scatter!(prof.log_R, LilGuys.middle.(prof.log_Sigma), label="recovered")
	errorbars!(prof.log_R, LilGuys.middle.(prof.log_Sigma), LilGuys.error_interval.(prof.log_Sigma))


	lines!(prof_nospace.log_R, prof_nospace.log_Sigma, label="true")

	lines!(prof_all.log_R, prof_all.log_Sigma, label="all")

	x = LinRange(-2, 2, 1000)
	R = 10 .^ x

	y = @. log10(Σ_bkd_model.(R) .* N_stars .* (1-f_sat) .* f_contam)
	lines!(x, y)
	text!(1, y[1], text="nospace bg", color=COLORS[4])

	limits!(-2, 2, -4, 4)
	axislegend(position=:lb)
	
	fig
end

# ╔═╡ 980992d4-93e2-4df3-81ca-bea1223d07b3
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(-2, 2, 1000)
	R = 10 .^ x
	ylims!(-5, 4)

	y = @. log10(Σ_memb.(R) .* f_sat .* N_stars)
	lines!(x, y, label="true initial density")


	y = @. log10((Σ_memb_model.(R) ./ Σ_memb_model(0)) .^ 0.1 .* f_sat .* N_stars .* Σ_memb.(R) )

	lines!(x, y, label="assumed density profile")
	
	
	y = @. log10(Σ_memb_model.(R) .* f_sat .* N_stars .* f_contam)

	lines!(x, y, label="assumed density profile contam")


	y = @. log10(Σ_bkd_model.(R) .* N_stars .* (1-f_sat))
	lines!(x, y)
	text!(1, y[1], text="total bg", color=COLORS[3])


	y = @. log10(Σ_bkd_model.(R) .* N_stars .* (1-f_sat) .* f_contam)
	lines!(x, y)
	text!(1, y[1], text="nospace bg", color=COLORS[4])


	scatter!(prof.log_R, LilGuys.middle.(prof.log_Sigma), label="recovered")
	errorbars!(prof.log_R, LilGuys.middle.(prof.log_Sigma), LilGuys.error_interval.(prof.log_Sigma))

	#axislegend(position=:lb)


	fig
end

# ╔═╡ 839a6688-9078-44a2-a314-15fd1db999cf
let
	fig = Figure()
	ax = Axis(fig[1,1])

	x = LinRange(-2, 2, 1000)
	R = 10 .^ x
	ylims!(0,1)

	y = @. (Σ_memb.(R) .* f_sat .* N_stars)
	lines!(x, y, label="true initial density")


	y = @. (Σ_memb_model.(R) .* f_sat .* N_stars )

	lines!(x, y, label="assumed density profile")


	y = @. (Σ_bkd_model.(R) .* N_stars .* (1-f_sat))
	lines!(x, y)
	text!(1, y[1], text="total bg", color=COLORS[3])


	y = @. (Σ_bkd_model.(R) .* N_stars .* (1-f_sat) .* f_contam)
	lines!(x, y)
	text!(1, y[1], text="nospace bg", color=COLORS[4])


	scatter!(prof.log_R, 10 .^ LilGuys.middle.(prof.log_Sigma), label="recovered")
	errorbars!(prof.log_R, 10 .^ LilGuys.middle.(prof.log_Sigma), LilGuys.error_interval.(prof.log_Sigma))

	#axislegend(position=:lb)


	fig
end

# ╔═╡ Cell order:
# ╠═f1b3eef7-7aae-4622-ba01-9470a2867b39
# ╠═446bb406-8ec0-4dba-822b-187c0b97d73b
# ╠═cc93ac36-795d-443e-998e-968a2d01f858
# ╠═a235a388-a42c-494f-8ef9-26554c5238b8
# ╠═82ff3270-fbe2-45df-943c-946f1292075a
# ╠═688e1152-81db-4a93-a314-1d902c87f054
# ╠═365909af-8fa6-49d1-8edb-3c7b646792e3
# ╠═34b0555d-605b-44e8-a467-8fb0962a54c4
# ╠═8e864e73-1c7f-4139-b921-db18fa004120
# ╠═d51390ea-f0cb-43f7-900a-539a901502f2
# ╠═0feed8ce-18c4-4213-8519-5749a4b16552
# ╠═ce09dc9f-0cea-4121-b0da-60a808399b2f
# ╠═1368181a-97d7-11f0-3b4f-b90907cb3133
# ╠═19112874-f462-4cfb-bbef-ca2a43e72fd1
# ╠═86671511-ac67-4837-bc29-1f4c8597519f
# ╠═a6af0647-c18e-4157-bf71-ed63186fe546
# ╠═2ee203cd-1370-441f-96de-bf086ca48bc6
# ╟─b719c6db-ec7c-45c3-a1d1-da9778c47477
# ╠═7288a3ad-375a-4930-a98b-2d225d8ed709
# ╠═99b70a40-b66b-4f24-b66e-f8e09aeca0f1
# ╠═7d78cfc8-18b9-416b-9676-1c86c9f84841
# ╠═b650ec10-1101-4032-8f42-34fb2b9681d2
# ╠═6ae2719c-f29a-4e6a-875e-90eed9c3a998
# ╠═5ed34292-2f95-4730-abcd-cb60a5509982
# ╠═db089c17-a106-4b61-9cff-8c6c04971355
# ╠═f2b73277-eccb-4e6c-8ec4-a76533e9fcdb
# ╠═7e827a7e-99fe-452d-be3e-24ba1a741887
# ╠═0afcd462-ab26-4b0c-a27a-3a838b0225e3
# ╠═d72a5173-5b76-4f30-8578-e615a5efecff
# ╠═9b918b65-0bf5-4564-8e09-c797adeb0b05
# ╠═f73aa823-08a9-4e21-bede-a9095647cd5a
# ╠═abfce12c-a3b2-4699-89db-734c87f1dd85
# ╠═59466065-3e1b-4cf2-9649-c48a10dfed05
# ╠═4fa3e12a-3407-4fcb-8452-b9cc46182d7f
# ╠═8021eb9e-9064-47d8-a108-afafa5c60ac4
# ╟─3f327a21-2532-4768-8fbb-4a68f2949a5f
# ╠═03a7b639-71a9-4c7f-9d97-57b07abbb3ae
# ╠═6aae7d18-eb3b-4ff8-96a1-06ef7a2ebfd8
# ╠═12711fe9-6518-455d-93e2-049d6cda0566
# ╠═bf609cbf-4755-4f43-9c8d-9f3601a29cf4
# ╠═5fc1c365-2634-4781-8729-7c4c774fb8cc
# ╠═bbb6cf08-48e4-4aa6-94a7-853c87d71b50
# ╠═3064eeaf-df33-40f3-8f61-2aa3405a7878
# ╠═4eace3ff-20fd-4a73-8e5e-8b9f9f5a34f5
# ╠═2675f07f-4a2b-46f8-a4a1-0fb994b6ef68
# ╠═a027b51c-aa7f-4158-80e3-c79f2c866e3d
# ╟─d6f5d1cb-2b0b-427e-b0f9-e2ef874a5eba
# ╠═10a4e8f4-8df5-49c3-b0f0-86b751936925
# ╠═4aed5645-e6a7-4935-8810-20e292149e77
# ╠═5c40c51f-690d-4d20-80ce-9446c220015b
# ╠═ed56f4f0-3946-4434-a613-59ef5efa974a
# ╠═1b8adcf9-2c19-4d5f-bb3c-5f6e77368a64
# ╠═61d5f5c3-daf7-4bd8-ac88-69a5ab453ab1
# ╠═e6b64d58-b9bb-4523-b38c-48625d7dabcb
# ╠═980992d4-93e2-4df3-81ca-bea1223d07b3
# ╠═839a6688-9078-44a2-a314-15fd1db999cf
