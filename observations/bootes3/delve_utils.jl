import LilGuys
using CSV, DataFrames
import KernelDensity
import Interpolations
import StatsBase: mean


struct DensityMap
    x::Vector{Float64}
    y::Vector{Float64}
    density::Matrix{Float64}
end

"""
	Retrieves the header of a Padova (PARSEC) formatted isochrone file.
Assumes that the header is the last row in the block of rows begining with 
`#` at the start of the file, with the columns seperated by spaces.
"""
function read_padova_header(filename)
	header = ""
	header_max = 30
	open(filename, "r") do f
		lastline = readline(f)
		for i in 1:header_max
			line = readline(f)
			if (line[1] != '#') && (lastline[1] == '#')
				header = lastline[3:end] # ignore comment character
				break
			end
			lastline = line
		end
	end
	header_vec = string.(split(header, r"\s+"))
end

function read_mist_file(filename)
	header = read_padova_header(filename)

	df = CSV.read(filename, DataFrame,
					comment="#", ignorerepeated=true, delim = ' ', header=header)
	return df
end


function get_isochrone(isochrones, age::Real, M_H::Real; stage_max::Integer=5)

	if age ∉ keys(isochrones)
		throw(KeyError("$age not in available ages"))
	end
	
	isos = isochrones[age]
	if (M_H < minimum(isos.MH)) || (M_H > maximum(isos.MH))
		throw(DomainError(M_H, "metallicity out of isochrone range: $(extrema(isos.MH))"))
	end

	M_Hs = unique(isos.MH)
	M_H_adopted = M_Hs[argmin(abs.(M_H .- M_Hs))]

	filt = isapprox.(isos.MH, M_H_adopted)
	filt .&= isos.label .<= stage_max # only keep through HB

	return isos[filt, :]
end



"""
    kroupa_imf(m) -> Float64

Kroupa (2001) broken power-law IMF. Returns dN/dm ∝ m^{-α}.
"""
function kroupa_imf(m::Real)
    m <= 0 && return 0.0
    if m < 0.08
        return m^(-0.3)
    elseif m < 0.50
        return 0.08^(1.3 - 0.3) * m^(-1.3)   # continuous at 0.08
    else
        return 0.08^(1.3 - 0.3) * 0.5^(2.3 - 1.3) * m^(-2.3)
    end
end

"""
    chabrier_imf(m) -> Float64

Chabrier (2003) IMF: log-normal below 1 M☉, Salpeter power-law above.
"""
function chabrier_imf(m::Real)
    m <= 0 && return 0.0
    if m <= 1.0
        return exp(-((log10(m) - log10(0.22))^2) / (2 * 0.57^2)) / m
    else
        return chabrier_imf(1.0) * m^(-2.3)
    end
end

# ╔═╡ bb0600b7-cf8a-484b-9d4a-7dab2c963da6
err_model(x, params) = @. (params[1] + x*params[2] + x^2 * params[3])

function fit_color_err(stars)

    gr_err = @. sqrt(stars.magerr_psf_g^2 + stars.magerr_psf_r^2)

    return fit_color_err(stars.gmag, gr_err)
end

function fit_color_err(gmag, gr_err)
    popt, covt = LilGuys.curve_fit(err_model, gmag, log10.(gr_err), [0., 1., 1.])

    @info popt
    color_err(gmag) = 10 .^ err_model(gmag, popt)

    return color_err
end



function resample_iso(iso; N=1000_000, cols=[:gmag, :rmag])
	massrange = extrema(iso.Mini)

	interps = Dict(col => LilGuys.lerp(iso.Mini, iso[!, col]) for col in cols)
	masses = LinRange(massrange..., N)
	df =  DataFrame(
		:Mini => masses,
	)

	for col in cols
		df[:, col] = interps[col].(masses)
	end
	df
end





"""
    build_cmd_filter(
        iso_color, iso_mag, iso_mass,
        color_uncertainty,
        distance_uncertainty;
        color_range  = (-0.5, 2.5),
        mag_range    = (15.0, 28.0),
        n_color      = 200,
        n_mag        = 300,
        imf          = kroupa_imf,
    ) -> (density, color_edges, mag_edges)

Build a matched-filter CMD density map from an isochrone, weighting
each point by the IMF and convolving along both the colour axis (photometric
uncertainty) and the magnitude axis (distance modulus uncertainty).

# Arguments
- `iso_color`            : Vector of isochrone colours (e.g. g-r)
- `iso_mag`              : Vector of isochrone magnitudes (same length)
- `iso_mass`             : Vector of stellar masses in solar units (same length)
- `color_uncertainty`    : Callable σ_c(mag) -> Float64, 1-σ colour error
- `distance_uncertainty` : Callable σ_μ(mag) -> Float64, 1-σ distance modulus
                           uncertainty in magnitudes. For a fixed σ_μ, pass
                           `_ -> σ_μ_value`. For a depth-dependent spread
                           (e.g. HB stars at different distances along the
                           line of sight), make this mag-dependent.

# Returns
- `density`     : Matrix{Float64} (n_color × n_mag), normalised to sum=1
- `color_edges` : Vector of colour bin edges (length n_color + 1)
- `mag_edges`   : Vector of magnitude bin edges (length n_mag + 1)
"""
function build_cmd_filter(
    iso_color::AbstractVector{<:Real},
    iso_mag::AbstractVector{<:Real},
    iso_mass::AbstractVector{<:Real},
    color_uncertainty,           # callable: σ_c(mag)   -> Float64
    distance_uncertainty;        # callable: σ_μ(mag)   -> Float64
    color_range::Tuple{Real,Real} = (-0.5, 1.5),
    mag_range::Tuple{Real,Real}   = (15.0, 23.0),
    n_color::Int = 200,
    n_mag::Int   = 300,
    imf = Utils.kroupa_imf,
    n_sigma_trunc=5,
)
    @assert length(iso_color) == length(iso_mag) == length(iso_mass) > 0

    # --- 1. Build output grid ---
    color_centres  = range(color_range[1], color_range[2]; length = n_color)
    mag_centres    = range(mag_range[1],   mag_range[2];   length = n_mag)
    Δc = step(color_centres)
    Δm = step(mag_centres)


    density = zeros(Float64, n_color, n_mag)

    # --- 2. Loop over isochrone points ---
    for k in 1:(length(iso_color)-1)
        c0 = (iso_color[k] + iso_color[k+1])/2
        m0 = (iso_mag[k] + iso_mag[k+1]) / 2
        ms = (iso_mass[k] + iso_mass[k+1]) / 2
		if (m0 > mag_range[2]) || (m0 < mag_range[1])
			continue
		end

        # IMF weight × mass interval
        Δmass = iso_mass[k+1] - iso_mass[k]
        w = imf(ms) * Δmass
        iszero(w) && continue

        # Uncertainty kernels at this point
        σ_c = max(color_uncertainty(m0),    Δc)   # colour σ, floored to 1 pixel
        σ_μ = max(distance_uncertainty(m0), Δm)   # distance modulus σ, floored

        # Precompute 1D kernel values to avoid redundant exp() calls
        kc = [LilGuys.gaussian(color_centres[i], c0, σ_c) for i in 1:n_color]

        # — magnitude kernel (length n_mag); only evaluate bins within ±4σ
        j_lo = max(1,    searchsortedlast(mag_centres, m0 - n_sigma_trunc*σ_μ))
        j_hi = min(n_mag, searchsortedlast(mag_centres, m0 + n_sigma_trunc*σ_μ))

        # --- 3. Outer product: colour kernel ⊗ magnitude kernel ---
        for j in j_lo:j_hi
            km = LilGuys.gaussian(mag_centres[j], m0, σ_μ)
            wkm = w * km
            iszero(wkm) && continue
            for i in 1:n_color
                density[i, j] += wkm * kc[i]
            end
        end
    end

    # --- 4. Normalise ---
    s = sum(density)
    s > 0 && (density ./= s)

    return DensityMap(collect(color_centres), collect(mag_centres), density) |> normalize_density
end



function normalize_density(dens_map::DensityMap)
    dx = LilGuys.gradient(dens_map.x) 
    dy = LilGuys.gradient(dens_map.y)'
    A = sum(dens_map.density .* dx .* dy)

    return DensityMap(dens_map.x, dens_map.y, dens_map.density ./ A)
end


function make_background_density(stars,     
        x, y, bandwidth = (0.05, 0.1)
)
    k = KernelDensity.kde((stars.gmag .- stars.rmag, stars.gmag), boundary=(extrema(x), extrema(y)), npoints=(length(x), length(y)), bandwidth=bandwidth)
    return DensityMap(k.x, k.y, k.density) |> normalize_density
end



function interpolate_density(dens_sat)
	dens = dens_sat.density 
	itp = Interpolations.interpolate((dens_sat.x, dens_sat.y), 
							   dens, Interpolations.Gridded(Interpolations.Linear()))	

	return Interpolations.extrapolate(itp, 0)
end


delve_g_err(gmag) = 10 ^ err_model(gmag, [-2.12775, -0.296242, 0.0147927])
delve_gr_err(gmag) = 10 ^ err_model(gmag, [-1.94422, -0.295023, 0.0146445])
gaia_G_err(gmag) =  10 ^ err_model(gmag, [3.1882, -0.927664, 0.0324749])
gaia_bp_rp_err(gmag) =  10 ^ err_model(gmag, [-3.5957, -0.137531, 0.0132695])



function create_kroupa_sampler(mass_max; mass_min=0.08)
	masses = LinRange(mass_min, mass_max, 10_000)
	imf =  kroupa_imf.(masses)
	N_cdf = cumsum(imf)
	x = N_cdf .- N_cdf[1]
	x ./= x[end]
	mass_func =  LilGuys.lerp(x, masses)

	function sampler()
		u = rand()
		return mass_func(u)
	end

end
