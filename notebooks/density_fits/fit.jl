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

# ╔═╡ 47b0d3e6-a79b-4f49-b847-e708e5d6aabf
md"""
 # setup
"""

# ╔═╡ acb9ae92-924b-4723-8bd7-d775595b24c3
COLORS = Arya.COLORS;

# ╔═╡ 95a463d3-fdec-4fa6-9f4f-6d4a485aebf1
F = Float64

# ╔═╡ d6cbd9cb-bfa7-46c1-970a-ab3fb3740c48
OptF = Union{F, Nothing}

# ╔═╡ e3f52a32-9a21-4c8f-8e5c-9ca3e7dbc331
begin
	Base.@kwdef struct DensityParams
		filename::String
		ra::Float64
		dec::Float64
		ecc::F
		rh::F
		PA::F
		dist::F
		dist_err::F
		PSAT_min::OptF = nothing
		ruwe_max::OptF = nothing
		g_min::OptF = nothing
		g_max::OptF = nothing
		max_ang_dist::OptF = nothing
		n_sigma_dist::OptF = nothing
		pmra::OptF = nothing
		pmdec::OptF = nothing
		dpm::OptF = nothing

		cmd_cut::Union{Array,Nothing} = nothing
	end

	function DensityParams(dict::Dict; kwargs...)
		d2 =  NamedTuple{Tuple(Symbol.(keys(dict)))}(values(dict))

		return DensityParams(; d2..., kwargs...)
	end
end

# ╔═╡ ff92927e-b078-45fd-9c13-1ce5a009d0bb
red = COLORS[6]

# ╔═╡ 8a551dbe-9112-48c2-be9a-8b688dc5a05c
md"""
# inputs
"""

# ╔═╡ 8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
name = "Scl.json"

# ╔═╡ 1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
begin 
	param_file = "params/$name"
	open(param_file, "r") do f
	    global params_json
	    params_json = JSON.parse(f)
	end

	params = DensityParams(params_json)
end

# ╔═╡ 4093a7d6-2f74-4c37-a4a8-270934ede924
md"""
# functions
"""

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
end

# ╔═╡ e4a6382e-21f6-4d34-8376-142ba859b18a
function Makie.convert_single_argument(y::Array{Measurement{T}}) where T
	return value.(y)
end

# ╔═╡ 695d532c-86d1-4b24-b7af-600a8ca29687
function plot_all_tangent!(ax, all_stars; scale=1, kwargs...)
    x = scale*all_stars.xi
    y = scale*all_stars.eta 
    return scatter!(ax, x, y; kwargs...)
end

# ╔═╡ 07235d51-10e1-4408-a4d1-cd2079fadb75
function plot_all_tangent(all_stars; scale=1, units="degrees", r_max=nothing, kwargs...)
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

# ╔═╡ 32fd9b79-1a8b-4a69-9115-9065dd61ced2
function load_fits(filename)
	f = FITS(filename)
	all_stars = DataFrame(f[2])
	close(f)
	return all_stars
end

# ╔═╡ c9fe69bf-52d5-4fec-9667-2c29d31230df
function add_r_ell!(stars, ecc, PA)
	b = sqrt(1 - ecc)
	a = 1/b
	
	println("a, b / r_h = $a, $b")
	r_ell = lguys.calc_r_ell(stars.xi, stars.eta, a, b, PA-90)
	stars[:, "r_ell"] = r_ell;
end

# ╔═╡ 0903be18-b239-4f38-98f2-2b170fc10c5a
function add_xi_eta!(stars, ra0, dec0)
	xi, eta = lguys.to_tangent(stars.ra, stars.dec, ra0, dec0)
	
	stars[:, "xi"] = xi
	stars[:, "eta"] = eta
	stars
end

# ╔═╡ 2ce9bedf-1ecb-4773-af65-1feff0f42f76
function min_filter(x, attr, cut)
	return x[:, attr] .> cut
end

# ╔═╡ ba162d65-2ee8-4390-9450-3975c649a05b
max_filter(x, attr, cut) = x[:, attr] .< cut

# ╔═╡ 2d8b2740-f9f2-4cca-823f-4b6491c31fe4
function apply_filter(df, func, params...)
	if any(params .== nothing)
		filt = trues(size(df, 1))
	else
		filt = func(df, params...)
	end
	println("filter cuts $func \t", sum(map(!, filt)))
	return filt
end

# ╔═╡ 60444f7f-71ee-4886-b02b-66bdc5324f99
md"""
# Filtering
"""

# ╔═╡ 029bb274-df79-4f8a-bdd2-13f293e38279
function cmd_filter(all_stars, cmd_cut)
	cmd_cut_m = reshape(cmd_cut, 2, :)
	filt_cmd = is_point_in_polygon.(zip(all_stars.bp_rp, all_stars.phot_g_mean_mag), [cmd_cut_m])
end

# ╔═╡ 4467fa31-7171-404e-a225-e3c27afa0f7d
function ang_dist_filter(all_stars, ra0, dec0, max_ang_dist)
filt_ang_dist = @. (
   	max_ang_dist ^2
    > (all_stars.ra - ra0)^2 * cosd(dec0)^2 
    + (all_stars.dec - dec0)^2
    )
end

# ╔═╡ d36ca6c6-f9bc-47c4-b728-e6badaf866b3
function pm_filter(all_stars, pmra, pmdec, dpm)
	filt_pm = @. (
	    dpm^2 
	    > (pmra - all_stars.pmra)^2 
	    + (pmdec - all_stars.pmdec)^2 
	    )
end

# ╔═╡ d86bb8bb-8d7a-4a22-964c-f3c8d90cc9f0
function psat_filter(all_stars, psat_min)
    println("number nan PSAT         \t", sum(isnan.(all_stars.PSAT)))
    println("number exactly zero PSAT\t", sum(all_stars.PSAT .== 0))
    println("number > zero           \t", sum(all_stars.PSAT .> 0))
    println("number == 1             \t", sum(all_stars.PSAT .== 1))

    println("total                   \t", length(all_stars.PSAT))
	return all_stars.PSAT .> psat_min
end

# ╔═╡ 9ee9ad05-41d2-4e26-8252-1ae322947bb1
function parallax_filter(all_stars, dist, dist_err, n_sigma_dist)
	parallax = 1/dist
	parallax_err = 1/dist * dist_err / dist_err

	sigma = @. sqrt(all_stars.parallax_error^2 + parallax_err^2)
	
	filt_parallax = @. (
    abs(all_stars.parallax - parallax) <  sigma * n_sigma_dist
    )
end

# ╔═╡ 914bdd71-da5c-4bf8-894d-12f64df5ca02
function select_members(all_stars, params)
	filt = apply_filter(all_stars, psat_filter, params.PSAT_min)

	filt .&= apply_filter(all_stars, min_filter, :phot_g_mean_mag, params.g_min)
	filt .&= apply_filter(all_stars, max_filter, :phot_g_mean_mag, params.g_max)
	filt .&= apply_filter(all_stars, max_filter, :ruwe, params.ruwe_max)
	filt .&= apply_filter(all_stars, cmd_filter, params.cmd_cut)
	filt .&= apply_filter(all_stars, ang_dist_filter, params.ra, params.dec, params.max_ang_dist)
	filt .&= apply_filter(all_stars, parallax_filter, params.dist, params.dist_err, params.n_sigma_dist)

	filt .&= apply_filter(all_stars, pm_filter, params.pmra, params.pmdec, params.dpm)



	println(sum(filt), " stars remaining")
	return all_stars[filt, :]
	
end

# ╔═╡ 753a7958-e12a-486e-b65b-7b6a8002a400
function load_and_filter(params)
	all_stars = load_fits(params.filename)
	add_xi_eta!(all_stars, params.ra, params.dec)
	add_r_ell!(all_stars, params.ecc, params.PA)

	members = select_members(all_stars, params)
	return all_stars, members
end

# ╔═╡ 44a44f97-9115-4610-9706-33acf065d0e7
all_stars, members = load_and_filter(params)

# ╔═╡ 0c498087-0184-4da2-a079-e972dd987712
md"""
The next three plots compare how different the r_ell and xi and eta calculated here are from what is (presumably) in the given catalogue.
"""

# ╔═╡ 52bb6b36-736a-45a8-b1e1-7f174b366ec8
let
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

# ╔═╡ f890216a-2e4e-4f44-92ee-ded0eaa17a68
params.rh

# ╔═╡ efc003db-c980-40ba-822f-23220f7e852e
md"""
# Membership plots
"""

# ╔═╡ d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
let 
	fig, ax, p = plot_all_tangent(all_stars, markersize=2,
        color=(:grey, 0.2))
	plot_all_tangent!(ax, members, markersize=2, color=red)
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	fig
end

# ╔═╡ 0010dffc-9717-4747-b7c2-2e396097399b
let 
	fig, ax, p = scatter(members.ra, members.dec, markersize=2,
        color=(:grey, 0.2))
	ax.aspect = 1 / cosd(params.dec)
	
	ax.xgridvisible = false
	ax.ygridvisible = false
	ax.xlabel="ra / degrees"
	ax.ylabel = "dec / degrees"
	fig
end

# ╔═╡ d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
let 
	fig, ax, p = plot_all_tangent(all_stars, markersize=5,
        color=(:grey, 0.2), r_max=30, scale=60, units="arcmin")
	plot_all_tangent!(ax, members, scale=60, markersize=5, color=red)
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
	    color=(red, 1), markersize=1)
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
	    color=(red, 1), markersize=1)
	
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
	    color=(red, 1), markersize=1)
	fig
end

# ╔═╡ 39f615b4-b62c-493e-8d11-be45910d79a8
md"""
# density calculation
"""

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
    
    return (log_r=log_r, log_r_bins=log_r_bin, counts=counts, Σ=Σ, Σ_m=Σ_m, Γ=Γ, Γ_max=Γ_max, M_in=M_in, N=length(rs))
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

# ╔═╡ 2c33af35-4217-4534-9e30-1a697b32892f


# ╔═╡ 5c117b2f-a32c-4afd-9c66-943ab4634e71
dist = params.dist ± params.dist_err

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

# ╔═╡ a7209445-84e9-435d-9144-90f8eb5e70cb
md"""
# Membership selection effects
"""

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

# ╔═╡ 45422d53-317c-4824-a41a-4a80b1fbd102
let 
	fig = Figure(size=(900, 500))
	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		xlabel=log_r_label,
	)

	g_cuts = [21, 20.5, 20, 19.5, 19, 10]
	labels = []

	ls = []
	for i in 1:length(g_cuts) - 1
		g_l = g_cuts[i+1]
		g_h = g_cuts[i]
		filt = g_cuts[i+1] .< members.phot_g_mean_mag .<= g_cuts[i]
		push!(labels, "$g_l, $g_h")
		memb = members[filt, :]
		println(size(memb, 1))
		r = memb.r_ell * 60
		obs = calc_properties(r)

		y = log10.(obs.Σ)
		f2 = isfinite.(y)
		y = y[f2]
		x = obs.log_r[f2]
		l = lines!(ax, x, value.(y), color=i, colorrange=(1, length(g_cuts) - 1))
		push!(ls, l)
	end

	Legend(fig[1,2], ls, labels, "G magnitude")

	fig
end

# ╔═╡ 13fb3ebc-50c0-43aa-88e9-1a7543e4e202
hist(members.phot_g_mean_mag)

# ╔═╡ 80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
size(members)

# ╔═╡ c0b3c3f6-0450-4242-9e13-41f9af17e562
let 
	fig = Figure(resolution=(900,500))
	p_cuts = [nothing, 0.01, 0.05, 0.1, 0.2, 0.5, 0.99]

	Nc = length(p_cuts) 

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(fraction/arcmin^2)}",
		xlabel=log_r_label,
		limits=(nothing, (-8, -1))
		
	)

	labels = string.(p_cuts)

	ls = []
	for i in 1:Nc
		params = DensityParams(params_json, PSAT_min=p_cuts[i])
		_, memb = load_and_filter(params)
		r = memb.r_ell * 60
		obs = calc_properties(r)
		l = lines!(ax, obs.log_r, log10.(value.(obs.Σ)), color=i, colorrange=(1, Nc))
		push!(ls, l)
	end

	Legend(fig[1,2], ls, labels, "probability cut")

	fig
end

# ╔═╡ 49ae0572-5d6b-4935-bc95-0a845bb3df2f
md"""
# Background density
"""

# ╔═╡ eeecad14-ec74-4559-9765-5648f0b3d74e
cmd_cut_umi =  [1.63,20.99, 1.26,19.88, 1.32,18.24, 1.69,16.24, 1.50,15.99, 1.19,17.31, 0.85,18.83, 0.67,19.29, -0.07,19.67, -0.18,20.08, 0.26,20.31, 0.58,20.11, 0.30,20.86]

# ╔═╡ d7984df8-84b1-41ff-b19b-dd17b1772d4a
r_max = maximum(sqrt.(all_stars.xi .^ 2 + all_stars.eta .^ 2))

# ╔═╡ 4fb45cd6-673e-48a6-a5a7-374abc7bf4a9
obs

# ╔═╡ c6362b4a-a2e8-4d3b-be03-79fb12b84b36
function scatter_dens!(obs, norm=:count; kwargs...)
	y = obs.Σ

	if norm == :count
		y *= obs.N
	end

	y = log10.(y)

	x = value.(obs.log_r)
	y_err = err.(y)
	y = value.(y)

	scatter!(x, y; kwargs...)
	errorbars!(x, y, y_err; kwargs...)
	
end
	

# ╔═╡ e033e344-737e-46e8-ab85-5fe33d191f41
"""
A simple density calculation 
"""
function calc_offset_density(dra, ddec, r_cut; n_sigma_dist=3, dpm=1, cmd_cut=cmd_cut_umi)

	
	params = DensityParams(params_json, 
		ra=params_json["ra"] + dra, dec=params_json["dec"] + ddec, PSAT_min=nothing, max_ang_dist=r_cut, ecc=0,
		dpm=dpm,
		cmd_cut=cmd_cut, n_sigma_dist=n_sigma_dist
	)

	_, memb = load_and_filter(params)
	r = memb.r_ell * 60
	
	obs = calc_properties(r)
	return obs
end

# ╔═╡ b6eaa6be-4a23-4357-9ce8-40aa9f16d7f6
let 
	fig = Figure(size=(900,500))

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label
	)

	r_shift = 0.75 * r_max
	r_cut = r_max - r_shift

	
	dras = r_shift*[1, 0, -1, 0, 0]
	ddecs = r_shift*[0, 1, 0, -1, 0]

	Nc = length(dras) 

	Σs = Measurement[]
	
	for i in 1:Nc
		dra = dras[i]
		ddec =  ddecs[i]
		label = "$(round(dra, digits=1)), $(round(ddec, digits=1))"
		obs = calc_offset_density(dra, ddec, r_cut)
		scatter_dens!(obs, label=label)
		if abs(dra^2 + ddec^2) > 0 
			append!(Σs, obs.Σ * obs.N)
		end
	end

	Σ_m = sum(Σs) / length(Σs)
	println("log Sigma background = $(log10.(Σ_m))")
	global log_Σ_bg
	log_Σ_bg = log10(Σ_m)
	hlines!(value.(log_Σ_bg))
	hspan!(value.(log_Σ_bg) .- err.(log_Σ_bg), value.(log_Σ_bg) .+ err.(log_Σ_bg), alpha=0.1)
	
	Legend(fig[1,2], ax, "position offset \n(degrees)", merge=true)

	fig
end

# ╔═╡ f832459e-edcb-48b4-ba3c-1d75a23f51e0
let 
	fig = Figure(size=(900,500))

	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label
	)

	r_shift = 0.75 * r_max
	r_cut = r_max - r_shift

	
	dras = r_shift*[1, 0, -1, 0, 0]
	ddecs = r_shift*[0, 1, 0, -1, 0]

	Nc = length(dras) 

	Σs = Measurement[]
	
	for i in 1:Nc
		dra = dras[i]
		ddec =  ddecs[i]
		label = "$(round(dra, digits=1)), $(round(ddec, digits=1))"
		obs = calc_offset_density(dra, ddec, r_cut, dpm=nothing, cmd_cut=nothing, n_sigma_dist=nothing)
		scatter_dens!(obs, label=label)
		if abs(dra^2 + ddec^2) > 0 
			append!(Σs, obs.Σ * obs.N)
		end
	end

	Σ_m = sum(Σs) / length(Σs)
	println("log Sigma background = $(log10.(Σ_m))")
	global log_Σ_bg2
	log_Σ_bg2 = log10(Σ_m)
	hlines!(value.(log_Σ_bg2))
	hspan!(value.(log_Σ_bg2) .- err.(log_Σ_bg2), value.(log_Σ_bg2) .+ err.(log_Σ_bg2), alpha=0.1)
	
	Legend(fig[1,2], ax, "position offset \n(degrees)", merge=true)

	fig
end

# ╔═╡ 48e41a6f-775d-4d9f-850d-df9bd20dcf09
let 
	fig = Figure(size=(900,500))


	ax = Axis(fig[1, 1],
		ylabel=L"\log \Sigma\ / \textrm{(number/arcmin^2)}",
		xlabel=log_r_label,
		
	)
	
	params = DensityParams(params_json)
	_, memb = load_and_filter(params)
	r = memb.r_ell * 60
	obs = calc_properties(r)
	scatter_dens!(obs, label="fiducial")

	obs = calc_offset_density(0, 0, 2)
	scatter_dens!(obs, label="simple")

	
	params = DensityParams(params_json, PSAT_min=nothing, ecc=0, )
	_, memb = load_and_filter(params)
	r = memb.r_ell * 60
	obs = calc_properties(r)
	scatter_dens!(obs, label="all")


	hlines!(value.(log_Σ_bg), label="background")
	hspan!(value.(log_Σ_bg) .- err.(log_Σ_bg), value.(log_Σ_bg) .+ err.(log_Σ_bg), alpha=0.1)

	hlines!(value.(log_Σ_bg2),  label="all backbround")
	hspan!(value.(log_Σ_bg2) .- err.(log_Σ_bg2), value.(log_Σ_bg2) .+ err.(log_Σ_bg2), alpha=0.1)
	
	Legend(fig[1, 2], ax, merge=true)
	fig
end

# ╔═╡ 44803049-4cc5-4a23-990a-322934ccb076
params_json

# ╔═╡ 45ce7a5d-75d4-4c7f-8233-5b2f7dde3a95
params

# ╔═╡ 28d71909-e3e7-4e32-9183-4da862a58eb7
begin 
	fits_name = name[begin:end-5] * "_fit.fits"
	f = FITS(fits_name, "w")
	obs_df = Dict{String, Any}()
	obs_units = Dict{String, Any}()

	
	obs_df["log_r"] = value.(obs.log_r)
	obs_df["log_r_err"] = err.(obs.log_r)
	obs_df["surface_dens"] = value.(obs.Σ)
	obs_df["surface_dens_err"] = err.(obs.Σ)
	obs_df["log_slope"] = value.(obs.Γ)
	obs_df["log_slope_err"] = err.(obs.Γ)

	obs_df["scale_radius_2d"] = value.(popt)
	obs_df["scale_radius_2d_err"] = err.(popt)


	
	obs_units["log_r"] = "log arcmin"
	obs_units["log_r_err"] = "log arcmin"
	obs_units["scale_radius_2d"] = "log arcmin"
	obs_units["scale_radius_2d_err"] = "log arcmin"
	
	write(f, obs_df, units=obs_units, hdutype=ASCIITableHDU)
	close(f)
end

# ╔═╡ Cell order:
# ╟─47b0d3e6-a79b-4f49-b847-e708e5d6aabf
# ╠═d5bec398-03e3-11ef-0930-f3bd4f3c64fd
# ╠═acb9ae92-924b-4723-8bd7-d775595b24c3
# ╠═95a463d3-fdec-4fa6-9f4f-6d4a485aebf1
# ╠═d6cbd9cb-bfa7-46c1-970a-ab3fb3740c48
# ╠═e3f52a32-9a21-4c8f-8e5c-9ca3e7dbc331
# ╠═ff92927e-b078-45fd-9c13-1ce5a009d0bb
# ╠═8a551dbe-9112-48c2-be9a-8b688dc5a05c
# ╠═8b2b3cec-baf7-4584-81bd-fa0a4fe2a4ac
# ╠═1514203c-8c64-49f2-bd2b-9b38e7e3e6ba
# ╠═4093a7d6-2f74-4c37-a4a8-270934ede924
# ╠═a0683e1e-5210-40fd-8841-1a6a315d3efe
# ╠═e4a6382e-21f6-4d34-8376-142ba859b18a
# ╠═4a9cf94b-8714-4320-9f8b-a480b0741ba5
# ╠═07235d51-10e1-4408-a4d1-cd2079fadb75
# ╠═695d532c-86d1-4b24-b7af-600a8ca29687
# ╠═32fd9b79-1a8b-4a69-9115-9065dd61ced2
# ╠═c9fe69bf-52d5-4fec-9667-2c29d31230df
# ╠═0903be18-b239-4f38-98f2-2b170fc10c5a
# ╠═44a44f97-9115-4610-9706-33acf065d0e7
# ╠═2ce9bedf-1ecb-4773-af65-1feff0f42f76
# ╠═ba162d65-2ee8-4390-9450-3975c649a05b
# ╠═2d8b2740-f9f2-4cca-823f-4b6491c31fe4
# ╠═753a7958-e12a-486e-b65b-7b6a8002a400
# ╟─60444f7f-71ee-4886-b02b-66bdc5324f99
# ╠═914bdd71-da5c-4bf8-894d-12f64df5ca02
# ╠═029bb274-df79-4f8a-bdd2-13f293e38279
# ╠═4467fa31-7171-404e-a225-e3c27afa0f7d
# ╠═d36ca6c6-f9bc-47c4-b728-e6badaf866b3
# ╠═d86bb8bb-8d7a-4a22-964c-f3c8d90cc9f0
# ╠═9ee9ad05-41d2-4e26-8252-1ae322947bb1
# ╟─0c498087-0184-4da2-a079-e972dd987712
# ╠═52bb6b36-736a-45a8-b1e1-7f174b366ec8
# ╠═f890216a-2e4e-4f44-92ee-ded0eaa17a68
# ╟─efc003db-c980-40ba-822f-23220f7e852e
# ╠═d7e51fb3-bfb2-4f19-963c-6a8eb497a88c
# ╠═0010dffc-9717-4747-b7c2-2e396097399b
# ╠═d0c4ad68-2bd3-44e7-9f42-45795fb7ecee
# ╠═b6424e6f-9b0d-4f29-b53d-0bd814a67139
# ╠═bffe11bd-4233-4a7c-9411-0dfb1ac79077
# ╠═0f002b56-8b8f-4025-8d7b-fb51423e8da0
# ╠═049ff11e-c04c-41d9-abf1-ec040b799649
# ╠═39f615b4-b62c-493e-8d11-be45910d79a8
# ╠═7c5397f3-40f3-49a4-acfb-27beadc5fa6b
# ╠═2d3308fa-ce87-4967-a934-d2821a42b8b7
# ╠═3d105351-9a5c-4410-a360-ec8736374909
# ╠═46bb8315-5bd7-4175-9da1-f8cad761b8ec
# ╠═7edaac34-e450-43a4-a022-90c38af9b5ca
# ╠═fb433202-4195-4dee-a2c6-f08656e64c2f
# ╠═830cbaab-b82f-4e8e-b975-1eb11662b11b
# ╠═4ceea6c3-0bf5-40c2-b49e-2691e73e003c
# ╠═58f7a246-f997-413b-abe4-73282abbc91c
# ╠═9b3288b4-3c17-4325-9e2c-94f96328f3c3
# ╠═1acc3b7d-2e9d-47ec-8afe-14876e57787c
# ╠═c0fa8744-4da3-470c-a2b8-f89ce431e1ed
# ╠═62388099-b80b-4272-9bde-0c7315b15c19
# ╠═e5d3d41b-e786-460a-a44f-51d1d98dcf11
# ╠═aca3780b-1d2f-4a44-a53a-754794c334b1
# ╠═2da8531e-eb37-4ec5-af5a-e0c2dfd8c4d1
# ╠═71d3b349-bf04-456a-b41c-a5acc58a7f73
# ╠═2c33af35-4217-4534-9e30-1a697b32892f
# ╠═5c117b2f-a32c-4afd-9c66-943ab4634e71
# ╠═0c87248a-cfa7-4a6d-af84-87f5543f68e2
# ╠═b20058a8-a9b6-49ff-b8ff-a2d45c76f645
# ╠═ac761271-6a6a-4df7-9476-1d44f02eb74d
# ╠═162dfe57-99e0-4c03-8934-2a59056f484f
# ╠═39cb37d9-f1d4-419e-9e19-c033bfba8556
# ╠═82d90218-f32e-4b72-a99a-bc2a264d7dce
# ╠═617a5128-e6a4-40a5-bc5e-45dba4eeaa58
# ╠═1ca6ea59-b129-460e-8925-2592332ba280
# ╟─a7209445-84e9-435d-9144-90f8eb5e70cb
# ╠═748430e6-520f-49c1-94e2-d1a0246d91b3
# ╠═45422d53-317c-4824-a41a-4a80b1fbd102
# ╠═13fb3ebc-50c0-43aa-88e9-1a7543e4e202
# ╠═80f2e2cf-c3b6-4931-b62f-4a2b9659fad5
# ╠═c0b3c3f6-0450-4242-9e13-41f9af17e562
# ╟─49ae0572-5d6b-4935-bc95-0a845bb3df2f
# ╠═eeecad14-ec74-4559-9765-5648f0b3d74e
# ╠═d7984df8-84b1-41ff-b19b-dd17b1772d4a
# ╠═4fb45cd6-673e-48a6-a5a7-374abc7bf4a9
# ╠═c6362b4a-a2e8-4d3b-be03-79fb12b84b36
# ╠═e033e344-737e-46e8-ab85-5fe33d191f41
# ╠═48e41a6f-775d-4d9f-850d-df9bd20dcf09
# ╠═b6eaa6be-4a23-4357-9ce8-40aa9f16d7f6
# ╠═f832459e-edcb-48b4-ba3c-1d75a23f51e0
# ╠═44803049-4cc5-4a23-990a-322934ccb076
# ╠═45ce7a5d-75d4-4c7f-8233-5b2f7dde3a95
# ╠═28d71909-e3e7-4e32-9183-4da862a58eb7
