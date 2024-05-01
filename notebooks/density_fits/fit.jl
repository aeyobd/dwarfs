### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ d5bec398-03e3-11ef-0930-f3bd4f3c64fd
begin 
	import Pkg; Pkg.activate()

	using FITSIO
	using DataFrames, CSV
	
	using GLMakie
	using Measurements
	#using KernelDensity
	
	import SciPy
	using QuadGK
	
	import LinearAlgebra: diag
	
	import LilGuys as lguys
	using Arya
	
	using JSON
end

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "UMi.json"

# ╔═╡ acb9ae92-924b-4723-8bd7-d775595b24c3
COLORS = Makie.wong_colors()

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
begin 
	param_file = "params/$name"
	params = Dict()
	open(param_file, "r") do f
	    global params
	    params = JSON.parse(f)
	end
end

# ╔═╡ 654d1641-6e14-46d7-a195-a050ff1c9980
outname = "results/$name"

# ╔═╡ 374e0815-b743-4033-a866-1b1a53b3b1b1
begin
	params["a"] = params["rh"] / 60
	params["b"] = (1-params["ecc"]) * params["a"]
end

# ╔═╡ 8a7b6d01-4d82-47d8-8cf8-56b990b62e21
fitsname = "data/" * params["filename"]

# ╔═╡ a0683e1e-5210-40fd-8841-1a6a315d3efe
function is_point_in_polygon(point, polygon)
    x, y = point
    inside = false
    N = size(polygon, 2)  # Number of vertices in the polygon
    j = N  # Start with the last vertex
    for i in 1:N
        xi, yi = polygon[1, i], polygon[2, i]
        xj, yj = polygon[1, j], polygon[2, j]
        
        # Check if point intersects with polygon edge
        intersect = ((yi > y) != (yj > y)) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
        if intersect
            inside = !inside
        end
        j = i
    end
    return inside
end

# ╔═╡ 4a9cf94b-8714-4320-9f8b-a480b0741ba5
begin 
	value(a::Measurement) = a.val
value(a) = a
err(a::Measurement) = a.err
err(a) = 0

function Makie.convert_single_argument(y::Array{Measurement{T}}) where T
    return value.(y)
end
end

# ╔═╡ 695d532c-86d1-4b24-b7af-600a8ca29687
function plot_all_tangent!(ax, all_stars; scale=1, kwargs...)
    x = scale*all_stars.xi
    y = scale*all_stars.eta 
    return scatter!(ax, x, y; kwargs...)
end

# ╔═╡ 07235d51-10e1-4408-a4d1-cd2079fadb75
function plot_all_tangent(all_stars; scale=1, units="degrees", r_max=nothing, kwargs...)
    #plot(xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}", aspect_ratio=1, xlim=(-10, 10), ylim=(-10, 10), dpi=400, fontfamily="Computer Modern", title="Gaia All")
    fig = Figure()

    x = scale*all_stars.xi
    y = scale*all_stars.eta
    
    if r_max === nothing
        r_max = max(maximum(abs.(x)), maximum(abs.(y)))
    end
    
    ax = Axis(fig[1,1], 
        xlabel=L"\xi / \textrm{%$units}", ylabel=L"\eta / \textrm{%$units}",
        aspect=1,
        limits=(-r_max, r_max, -r_max, r_max)
    )

    p = plot_all_tangent!(ax, all_stars; scale=scale, kwargs...) 

    return Makie.FigureAxisPlot(fig, ax, p)
end

# ╔═╡ 8eb9d6ee-8d10-4556-a2d2-093738a99869
begin
	f = FITS(fitsname)
	all_stars = DataFrame(f[2]);

	xi, eta = lguys.to_tangent(all_stars.ra, all_stars.dec, params["ra"], params["dec"],)
	r_ell = lguys.calc_r_ell(xi, eta, params["a"], params["b"], params["PA"]-90) * sqrt(params["a"]*params["b"])
	
	all_stars[:, "xi_0"] = copy(all_stars.xi)
	all_stars[:, "eta_0"] = copy(all_stars.eta)
	all_stars[:, "r_ell_0"] = copy(all_stars.r_ell)
	
	all_stars[:, "xi"] = xi
	all_stars[:, "eta"] = eta
	all_stars[:, "r_ell"] = r_ell;
end

# ╔═╡ 60444f7f-71ee-4886-b02b-66bdc5324f99
md"""
# Filtering
"""

# ╔═╡ 91836e7a-3812-4a8b-aa2a-81b2524c6f37
jax = params["filter_mode"] == "jax"

# ╔═╡ 1ecfa947-3376-4695-a042-24c1e46db12d
begin 
	params["parallax"] = 1 / params["dist"]
params["parallax_err"] = params["dist_err"] / params["dist"] * params["parallax"]
end

# ╔═╡ 52bb6b36-736a-45a8-b1e1-7f174b366ec8
let
	if params["filter_mode"] == "jax"
	    f = Figure()
	    ax = Axis(f[1, 1],
	        xlabel="probability", 
	        ylabel="count",
	        yscale=log10)
	
	    hist!(ax, all_stars.PSAT[all_stars.PSAT .>= 0], 
	        bins=20, label="j24")
	    #stephist!(b22.Pmemb, label="b22")
		f
	end
end

# ╔═╡ 028cd5e5-3475-4e0b-9bc7-c793ad11376d
if params["filter_mode"] == "jax"

    println("number nan PSAT         \t", sum(isnan.(all_stars.PSAT)))
    println("number exactly zero PSAT\t", sum(all_stars.PSAT .== 0))
    println("number > zero           \t", sum(all_stars.PSAT .> 0))
    println("number == 1             \t", sum(all_stars.PSAT .== 1))

    println("total                   \t", length(all_stars.PSAT))
    filt_jax = all_stars.PSAT .>  params["psat_min"]
 
end;

# ╔═╡ d1685148-72a2-436e-b896-07fb82ab3dcf
md"""
## Simple filtering
"""

# ╔═╡ 78e150e1-8855-490d-ab0c-14a0a6d49f59
if params["filter_mode"] == "simple"
    cmd_cut = reshape(params["cmd_cut"], 2, :)
	
	filt_cmd = is_point_in_polygon.(zip(all_stars.bp_rp, all_stars.phot_g_mean_mag), [cmd_cut])
	println(sum(filt_cmd))
end

# ╔═╡ 232ac5b5-eb7f-430a-9ffd-051d245a2362

filt_ang_dist = @. (
    params["max_ang_dist"]
    > (all_stars.ra - params["ra"])^2 * cosd(params["dec"])^2 
    + (all_stars.dec - params["dec"])^2
    )

# ╔═╡ c424cd66-c7c6-40c9-bd35-99dbbc75d1e4
filt_pm = @. (
    params["dpm"]^2 
    > (params["pmra"] - all_stars.pmra)^2 
    + (params["pmdec"] - all_stars.pmdec)^2 
    )

# ╔═╡ c3e72bdd-1cd0-4ec4-890c-f70a5f1f254d
filt_parallax = @. (
    abs(all_stars.parallax - params["parallax"]) 
    < params["n_sigma_dist"]
    * sqrt(all_stars.parallax_error^2 + params["parallax_err"]^2)
    )

# ╔═╡ a4be4772-fb1e-44a2-b18a-dba3c1934710
begin 
	filt_good = all_stars.ruwe .< params["ruwe_max"]
filt_good .&= all_stars.phot_g_mean_mag .< params["g_min"];
filt_good .&= all_stars.phot_g_mean_mag .> params["g_max"];
end

# ╔═╡ 35d56357-48a4-44e8-b099-c024de5c74ab
if params["filter_mode"] == "simple"
	filt_me = filt_parallax .& filt_pm .& filt_cmd .& filt_good .& filt_ang_dist;
end

# ╔═╡ 4becf3bb-c5b1-49e6-a6ba-d87160f80266
md"""
## Appluing filter
"""

# ╔═╡ 6e0ac8e1-b019-4cde-933b-f82a96b93a73
if params["filter_mode"] == "jax"
    filt = filt_jax
elseif params["filter_mode"] == "simple"
    filt = filt_me
else
end

# ╔═╡ 696d69ee-fedc-4b94-9026-45928891c8d0
members = all_stars[filt, :]

# ╔═╡ 0c498087-0184-4da2-a079-e972dd987712
md"""
The next three plots compare how different the r_ell and xi and eta calculated here are from what is (presumably) in the given catalogue.
"""

# ╔═╡ e8092e03-359f-4560-9c3e-130c94a1803b
if jax
	let
		fig, ax, p = hist(all_stars.xi .- all_stars.xi_0)
		ax.xticklabelrotation = -π/5
		ax.xlabel = L"\Delta \xi"
		fig
	end
end

# ╔═╡ 6e329265-614a-40af-998d-af9b6aa30f70
if jax
let
	fig, ax, p = hist(all_stars.eta .- all_stars.eta_0)
	ax.xticklabelrotation = -π/5
	ax.xlabel = L"\Delta \eta"
	fig
end
end

# ╔═╡ eeb06a82-d2ad-46d4-9d4d-832ada452fbe
if jax
let
	fig, ax, p = hist(all_stars.r_ell .- all_stars.r_ell_0 * sqrt(params["a"]*params["b"]))

	ax.xlabel = L"\Delta r_\textrm{ell}"
	fig
end
end

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
let 
	fig, ax, p = plot_all_tangent(all_stars, markersize=2,
        color=(:grey, 0.2))
	plot_all_tangent!(ax, members, markersize=2, color=(COLORS[6]))
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	fig
end

# ╔═╡ d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
let 
	fig, ax, p = plot_all_tangent(all_stars, markersize=5,
        color=(:grey, 0.2), r_max=30, scale=60, units="arcmin")
	plot_all_tangent!(ax, members, scale=60, markersize=5, color=COLORS[6])
	fig
end

# ╔═╡ b6424e6f-9b0d-4f29-b53d-0bd814a67139
let	
	fig = Figure()
	da = 60
	ax = Axis(fig[1, 1], 
	    xlabel=L"\xi / \textrm{arcmin} ", ylabel=L"\eta / \textrm{arcmin}",
	    aspect=1,
	    limits = (-da, da, -da, da))
	
	scatter!(ax, 
	    60*members.xi, 60*members.eta, 
	    color=members.phot_g_mean_mag, 
	    colormap=:greys,
	    markersize=5
	)
	
	
	fig
end

# ╔═╡ bffe11bd-4233-4a7c-9411-0dfb1ac79077
let
	fig = Figure()
	da = 15
	ax = Axis(fig[1, 1], limits=(-da, da, -da, da), aspect=1,
	    xlabel=L"\mu_{\alpha *} / \mathrm{ mas\, yr^{-1}}", 
	    ylabel=L"\mu_\delta / \mathrm{ mas\, yr^{-1}}")
	scatter!(ax, all_stars.pmra, all_stars.pmdec, 
	    color=(:grey, 0.2), markersize=1)
	
	scatter!(ax, members.pmra, members.pmdec, 
	    color=(:red, 1), markersize=1)
	fig
end

# ╔═╡ 0f002b56-8b8f-4025-8d7b-fb51423e8da0
let
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=1,
	    limits=(-0.5, 3, 10, 22), yreversed=true,
	    xlabel="bp - rp", 
	    ylabel="G",)
	scatter!(ax, all_stars.bp_rp, all_stars.phot_g_mean_mag, 
	    color=(:grey, 0.2), markersize=1)
	
	scatter!(ax, members.bp_rp, members.phot_g_mean_mag, 
	    color=(:red, 1), markersize=1)
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	
	fig
end

# ╔═╡ 049ff11e-c04c-41d9-abf1-ec040b799649
let
	fig = Figure()
	da = 15
	ax = Axis(fig[1, 1], aspect=1, limits=(-10, 10, 0, 4),
	    xlabel=L"\varpi / \mathrm{ mas }", 
	    ylabel=L"\delta\varpi / \mathrm{ mas }")
	scatter!(ax, all_stars.parallax, all_stars.parallax_error, 
	    color=(:grey, 0.2), markersize=1)
	
	scatter!(ax, members.parallax, members.parallax_error, 
	    color=(:red, 1), markersize=1)
	fig
end

# ╔═╡ 7c5397f3-40f3-49a4-acfb-27beadc5fa6b
function running_hist(xs, bw, normalize=false)
	N = length(xs)
	hist = zeros(N)

	x_sort = sort(xs)
	for i in 1:N
		x0 = x_sort[i]
		count = sum(x0 + bw .> x_sort .> x0 - bw)
		hist[i] = count
	end

	if normalize
		dx = lguys.gradient(x_sort)
		area = sum(hist .* dx)
		hist ./= area
	end
	return x_sort, hist
end

# ╔═╡ 2d3308fa-ce87-4967-a934-d2821a42b8b7
function norm_hist(xs, bw)
	bins = collect(minimum(xs):bw:(maximum(xs) + bw))
	x_h, y_h = lguys.calc_histogram(xs, bins)
    y_e = sqrt.(y_h)
    
	area = length(xs)

	return x_h, y_h ./ area, y_e ./ area
end

# ╔═╡ 3d105351-9a5c-4410-a360-ec8736374909
function calc_Σ(log_r, hist)
	r = 10 .^ log_r
	Σ = hist ./ (2π * log(10) * r .^ 2) # is equivalent because of derivative of log r
	return Σ
end

# ╔═╡ 46bb8315-5bd7-4175-9da1-f8cad761b8ec
function calc_Σ_mean(log_r, hist)
	r = 10 .^ lguys.midpoint(log_r)
	counts = cumsum(hist .* diff(log_r))
	Areas = @. π * r^2
	σ = counts ./ Areas
	return σ
end

# ╔═╡ 7edaac34-e450-43a4-a022-90c38af9b5ca
function calc_Γ(log_rs, ρs, step=1)
	dx = lguys.gradient(log_rs)
	dρ = lguys.gradient(log10.(ρs))

	return dρ ./ dx #lguys.gradient(log10.(ρs), log_rs)
end

# ╔═╡ fb433202-4195-4dee-a2c6-f08656e64c2f
function log_Σ_exp(r, A, r_s)
    return @. A + 1/log(10) * (-r / r_s)
end

# ╔═╡ 830cbaab-b82f-4e8e-b975-1eb11662b11b
function calc_properties(rs, bw=0.1)
    log_r_bin, counts, δ_counts = norm_hist(log10.(rs), 0.1)
    counts = counts .± δ_counts
    log_r = lguys.midpoint(log_r_bin)
    δ_log_r = diff(log_r_bin) ./ 2
    log_r = log_r .± δ_log_r
    r_bin = 10 .^ log_r_bin

    r = 10 .^ log_r

    As = π * diff((10 .^ log_r_bin) .^ 2)
    Σ = counts ./ As 
    # Σ_e = @. y_e / ys * Σ
    #Σ_m = calc_Σ_mean(xs, ys)

    M_in = cumsum(counts)
    A_in = @. π * (r_bin[2:end])^2
    Σ_m = M_in ./ A_in

    Γ = calc_Γ(log_r, Σ)
    Γ_max = @. 2*(1 - Σ / Σ_m)
    
    return (log_r=log_r, log_r_bins=log_r_bin, counts=counts, Σ=Σ, Σ_m=Σ_m, Γ=Γ, Γ_max=Γ_max, M_in=M_in)
end

# ╔═╡ 58f7a246-f997-413b-abe4-73282abbc91c
function predict_properties(Σ_model; N=10_000, log_r_min=-2, log_r_max=2)
    log_r_bins = LinRange(log_r_min, log_r_max, 1000)
    log_r = lguys.midpoint(log_r_bins)
    r = 10 .^  log_r
    r_bins = 10 .^ log_r_bins
    
    Σ = Σ_model.(r)
    Γ = calc_Γ(log_r, Σ)
    M_in = [quadgk(rrr->2π*rrr*Σ_model(rrr), 0, rr)[1] for rr in r]
    Σ_m = M_in ./ (π * r .^ 2)
    Γ_max = 2*(1 .- Σ ./ Σ_m)
    counts =  Σ .* (2π *  r .* diff(r_bins) )
    
    return (log_r=log_r, log_r_bins=log_r_bins, counts=counts, Σ=Σ, Σ_m=Σ_m, Γ=Γ, Γ_max=Γ_max, M_in=M_in)

end

# ╔═╡ 4ceea6c3-0bf5-40c2-b49e-2691e73e003c
function fit_profile(obs; r_max=Inf, N=10_000)
    r_val = [10 ^ r.val for r in obs.log_r]
    log_Σ = log10.(obs.Σ)
    filt = r_val .< r_max
    filt .&= map(x->isfinite(x), log_Σ)
    filt .&= @. !isnan(log_Σ)
    
    r_val = r_val[filt]
    log_Σ = log_Σ[filt]


    log_Σ_val = [s.val for s in log_Σ]
    log_Σ_e = [s.err for s in log_Σ]

    popt, covt = SciPy.optimize.curve_fit(log_Σ_exp, r_val, log_Σ_val, 
        sigma=log_Σ_e, p0=[1, 0.1])
    
    popt_p = popt .± sqrt.(diag(covt))
    println("log_Σ_0 = $(popt_p[1])")
    println("r_s = $(popt_p[2])")
    props = predict_properties(r->10 .^ log_Σ_exp(r, popt...), 
        N=N, log_r_min=obs.log_r_bins[1], log_r_max=obs.log_r_bins[end])
    
    log_Σ_pred = log_Σ_exp.(10 .^ value.(obs.log_r), popt...)
    log_Σ_res = log10.(obs.Σ) .- log_Σ_pred
    return popt_p, props, log_Σ_res
end

# ╔═╡ 3d1950d9-1fad-4002-9bd0-1cc517a7bbea
"""
    normalize_histogram(bins, counts)

Normalizes a histogram such that the area is equal to 1, return the bin midpoints
and the normalized (pdf) of the counts
"""
function normalize_histogram(bins, counts)
    A = sum(counts .* diff(bins))
    
    return lguys.midpoint(bins), counts ./ A
end

# ╔═╡ 9b3288b4-3c17-4325-9e2c-94f96328f3c3
function plot_rh!()
    vline!([log10.(r_h)], color="grey", s=:dash, z_order=1, label=L"r_h")
end

# ╔═╡ c0fa8744-4da3-470c-a2b8-f89ce431e1ed
log_r_label = L"\log r / \mathrm{arcmin}"

# ╔═╡ 1acc3b7d-2e9d-47ec-8afe-14876e57787c
function plot_Σ_fit_res(obs, pred, res)
    fig = Figure()
    ax = Axis(fig[1, 1], 
        ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}")
    N = length(obs.log_r)
    y = log10.(obs.Σ .* N)
    errorbars!(ax, value.(obs.log_r), value.(y), err.(y))
    scatter!(ax, value.(obs.log_r), value.(y))

    lines!(ax, pred.log_r, log10.(pred.Σ .* N), color=COLORS[2])
    
    ax2 = Axis(fig[2, 1],
        ylabel=L"\delta\log\Sigma", 
    	xlabel=log_r_label,
		limits = (nothing, (-1, 1))
	)

#     p2 = plot(ylabel=L"\Delta \log\Sigma", xlabel=log_r_label, ylim=(-2, 2))

	y = res
    scatter!(ax2, value.(obs.log_r), value.(y), err.(y), msw=1, msc=1, label="")
    errorbars!(ax2, value.(obs.log_r), value.(y), err.(y))


    hlines!(0, color=:black)
    
    rowsize!(fig.layout, 2, Relative(1/4))

    linkxaxes!(ax, ax2)
    hidexdecorations!(ax, grid=false)
#     return plot(p1, p2, layout=grid(2, 1, heights=(0.8, 0.2)), link=:x, bottom_margin=[-5Plots.mm 0Plots.mm])
    return fig
end

# ╔═╡ 62388099-b80b-4272-9bde-0c7315b15c19
r = 60 * members.r_ell # arcminutes

# ╔═╡ e5d3d41b-e786-460a-a44f-51d1d98dcf11
obs = calc_properties(r)

# ╔═╡ aca3780b-1d2f-4a44-a53a-754794c334b1
println(maximum(r) / 60)

# ╔═╡ 2da8531e-eb37-4ec5-af5a-e0c2dfd8c4d1
let
	fig, ax, p = hist(log10.(r))
	ax.xlabel = log_r_label
	ax.ylabel = "count"
	fig
end

# ╔═╡ 71d3b349-bf04-456a-b41c-a5acc58a7f73
popt, pred, res = fit_profile(obs)

# ╔═╡ 5c117b2f-a32c-4afd-9c66-943ab4634e71
dist = params["dist"] ± params["dist_err"]

# ╔═╡ b20058a8-a9b6-49ff-b8ff-a2d45c76f645
R_s_over_R_h = 1.6783

# ╔═╡ 0c87248a-cfa7-4a6d-af84-87f5543f68e2
popt[2]  * 60 / (206265) * dist * 1e3 * R_s_over_R_h # scale radius in pc

# ╔═╡ ac761271-6a6a-4df7-9476-1d44f02eb74d
popt[2] * 1.6783 # half light radius

# ╔═╡ 162dfe57-99e0-4c03-8934-2a59056f484f
plot_Σ_fit_res(obs, pred, res)

# ╔═╡ 39cb37d9-f1d4-419e-9e19-c033bfba8556
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=(-0.8, 1.7, -8, 5))
	scatter!(value.(obs.log_r), value.(obs.Γ),)
	errorbars!(value.(obs.log_r), value.(obs.Γ), err.(obs.Γ) )
	
	lines!(pred.log_r, pred.Γ, label="exponential", color=COLORS[2])
	
	ax.xlabel = log_r_label
	ax.ylabel = L"\Gamma = d\,\log \Sigma / d\,\log r"
	
	fig
end

# ╔═╡ 82d90218-f32e-4b72-a99a-bc2a264d7dce
theme(:colorcycle)

# ╔═╡ 617a5128-e6a4-40a5-bc5e-45dba4eeaa58
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	scatter!(10 .^ value.(obs.log_r), value.(obs.Γ),)
	errorbars!(10 .^ value.(obs.log_r), value.(obs.Γ), err.(obs.Γ) )
	
	lines!(10 .^ pred.log_r, pred.Γ, label="exponential", color=COLORS[2])
	
	ax.xlabel = "r / arcmin"
	ax.ylabel = L"\Gamma = d\,\log \Sigma / d\,\log r"
	
	fig
end

# ╔═╡ 1ca6ea59-b129-460e-8925-2592332ba280
let
	fig = Figure()
	ax = Axis(fig[1, 1], limits=(-0.8, 2.2, -1, 2.2))
	scatter!(value.(obs.log_r), value.(obs.Γ_max),)
	errorbars!(value.(obs.log_r), value.(obs.Γ_max), err.(obs.Γ_max) )
	
	lines!(pred.log_r, pred.Γ_max, label="exponential", color=COLORS[2])
	
	ax.xlabel = log_r_label
	ax.ylabel = L"\Gamma_\mathrm{max} = 2(1 - \Sigma / \bar{\Sigma})"
	
	fig
end

# ╔═╡ 748430e6-520f-49c1-94e2-d1a0246d91b3
begin 
	data = Dict{String, Any}()
	
	data["log_r"] = value.(obs.log_r)
	data["Sigma"] = value.(obs.Σ)
	data["Sigma_err"] = err.(obs.Σ)
	data["Gamma"] = value.(obs.Γ)
	data["Gamma_err"] = err.(obs.Γ)
	data["log_Sigma_0"] = popt[1].val
	data["log_Sigma_0_err"] = popt[1].err
	data["r_s"] = popt[2].val
	data["r_s_err"] = popt[2].err
	data
end

# ╔═╡ c7bb70d5-f2dc-45f8-af0a-67b981b10f88
open(outname, "w") do f
    s = json(data)
    write(f, s)
end

# ╔═╡ Cell order:
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═acb9ae92-924b-4723-8bd7-d775595b24c3
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╠═654d1641-6e14-46d7-a195-a050ff1c9980
# ╠═374e0815-b743-4033-a866-1b1a53b3b1b1
# ╠═8a7b6d01-4d82-47d8-8cf8-56b990b62e21
# ╠═a0683e1e-5210-40fd-8841-1a6a315d3efe
# ╠═4a9cf94b-8714-4320-9f8b-a480b0741ba5
# ╠═07235d51-10e1-4408-a4d1-cd2079fadb75
# ╠═695d532c-86d1-4b24-b7af-600a8ca29687
# ╠═8eb9d6ee-8d10-4556-a2d2-093738a99869
# ╠═60444f7f-71ee-4886-b02b-66bdc5324f99
# ╠═91836e7a-3812-4a8b-aa2a-81b2524c6f37
# ╠═1ecfa947-3376-4695-a042-24c1e46db12d
# ╠═52bb6b36-736a-45a8-b1e1-7f174b366ec8
# ╠═028cd5e5-3475-4e0b-9bc7-c793ad11376d
# ╠═d1685148-72a2-436e-b896-07fb82ab3dcf
# ╠═78e150e1-8855-490d-ab0c-14a0a6d49f59
# ╠═232ac5b5-eb7f-430a-9ffd-051d245a2362
# ╠═c424cd66-c7c6-40c9-bd35-99dbbc75d1e4
# ╠═c3e72bdd-1cd0-4ec4-890c-f70a5f1f254d
# ╠═a4be4772-fb1e-44a2-b18a-dba3c1934710
# ╠═35d56357-48a4-44e8-b099-c024de5c74ab
# ╠═4becf3bb-c5b1-49e6-a6ba-d87160f80266
# ╠═6e0ac8e1-b019-4cde-933b-f82a96b93a73
# ╠═696d69ee-fedc-4b94-9026-45928891c8d0
# ╟─0c498087-0184-4da2-a079-e972dd987712
# ╠═e8092e03-359f-4560-9c3e-130c94a1803b
# ╠═6e329265-614a-40af-998d-af9b6aa30f70
# ╠═eeb06a82-d2ad-46d4-9d4d-832ada452fbe
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
# ╠═d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
# ╠═b6424e6f-9b0d-4f29-b53d-0bd814a67139
# ╠═bffe11bd-4233-4a7c-9411-0dfb1ac79077
# ╠═0f002b56-8b8f-4025-8d7b-fb51423e8da0
# ╠═049ff11e-c04c-41d9-abf1-ec040b799649
# ╠═7c5397f3-40f3-49a4-acfb-27beadc5fa6b
# ╠═2d3308fa-ce87-4967-a934-d2821a42b8b7
# ╠═3d105351-9a5c-4410-a360-ec8736374909
# ╠═46bb8315-5bd7-4175-9da1-f8cad761b8ec
# ╠═7edaac34-e450-43a4-a022-90c38af9b5ca
# ╠═fb433202-4195-4dee-a2c6-f08656e64c2f
# ╠═830cbaab-b82f-4e8e-b975-1eb11662b11b
# ╠═4ceea6c3-0bf5-40c2-b49e-2691e73e003c
# ╠═58f7a246-f997-413b-abe4-73282abbc91c
# ╠═3d1950d9-1fad-4002-9bd0-1cc517a7bbea
# ╠═9b3288b4-3c17-4325-9e2c-94f96328f3c3
# ╠═1acc3b7d-2e9d-47ec-8afe-14876e57787c
# ╠═c0fa8744-4da3-470c-a2b8-f89ce431e1ed
# ╠═62388099-b80b-4272-9bde-0c7315b15c19
# ╠═e5d3d41b-e786-460a-a44f-51d1d98dcf11
# ╠═aca3780b-1d2f-4a44-a53a-754794c334b1
# ╠═2da8531e-eb37-4ec5-af5a-e0c2dfd8c4d1
# ╠═71d3b349-bf04-456a-b41c-a5acc58a7f73
# ╠═5c117b2f-a32c-4afd-9c66-943ab4634e71
# ╠═0c87248a-cfa7-4a6d-af84-87f5543f68e2
# ╠═b20058a8-a9b6-49ff-b8ff-a2d45c76f645
# ╠═ac761271-6a6a-4df7-9476-1d44f02eb74d
# ╠═162dfe57-99e0-4c03-8934-2a59056f484f
# ╠═39cb37d9-f1d4-419e-9e19-c033bfba8556
# ╠═82d90218-f32e-4b72-a99a-bc2a264d7dce
# ╠═617a5128-e6a4-40a5-bc5e-45dba4eeaa58
# ╠═1ca6ea59-b129-460e-8925-2592332ba280
# ╠═748430e6-520f-49c1-94e2-d1a0246d91b3
# ╠═c7bb70d5-f2dc-45f8-af0a-67b981b10f88
