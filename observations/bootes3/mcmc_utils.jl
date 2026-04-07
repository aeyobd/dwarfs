using LilGuys
using Turing
using DataFrames
import StatsBase: quantile


FVEC = AbstractVector{<:Real}
const LL_min = 1e-300
const PVALUE = 0.16
const α_exp = LilGuys.R_h(LilGuys.Exp2D())

prior_f_sat = Uniform(0, 0.5)
prior_R_h = LogNormal(log(30), 0.5)
prior_tangent = Normal(0, 15)

Base.@kwdef struct GaiaData
    source_id::Vector{Int}
    xi::Vector{Float64}
    eta::Vector{Float64}
    "Prior BG probability"
    L_bg::Vector{Float64}
    "Prior Satellite probability"
    L_sat::Vector{Float64}
end

@model function exp_ell_model(data::GaiaData;
        area_tot::Real,
        prior_tangent = prior_tangent, 
        prior_R_h = prior_R_h, 
        prior_f_sat = prior_f_sat
    )
	d_xi  ~ prior_tangent
	d_eta ~ prior_tangent
	R_h ~ prior_R_h
	f_sat ~ prior_f_sat

	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	
	R = LilGuys.calc_R_ell(data.xi .- d_xi, data.eta .- d_eta, 
                           ellipticity, position_angle)

    prof = LilGuys.Exp2D(R_s=R_h/α_exp, M=1)

	L_sat_space = @. LilGuys.surface_density.(prof, R) 
	L_bg_space = 1 / area_tot

    L = @. ((1-f_sat) * L_bg_space * data.L_bg
            + f_sat *  L_sat_space * data.L_sat)
    L = max.(L, LL_min)
    LL = sum(log.(L))
	
	Turing.@addlogprob!(LL)
end


@model function plummer_model(data::GaiaData;
        area_tot::Real,
        prior_tangent = prior_tangent, 
        prior_R_h = prior_R_h, 
        prior_f_sat = prior_f_sat
    )
	d_xi  ~ prior_tangent
	d_eta ~ prior_tangent
	R_h ~ prior_R_h
	f_sat ~ prior_f_sat

	R = @. sqrt((data.xi - d_xi)^2 + (data.eta - d_eta)^2)

	prof = LilGuys.Plummer(r_s=R_h, M=1)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R) 
	L_bg_space = 1 / area_tot
	
	
    L = @. ((1-f_sat) * L_bg_space * data.L_bg
            + f_sat *  L_sat_space * data.L_sat)
    L = max.(L, LL_min)
    LL = sum(log.(L))
	
	Turing.@addlogprob!(LL)
end


@model function plummer_ell_model(data::GaiaData;
        area_tot::Real,
        prior_tangent = prior_tangent, 
        prior_R_h = prior_R_h, 
        prior_f_sat = prior_f_sat
    )

	d_xi  ~ prior_tangent
	d_eta ~ prior_tangent
	R_h ~ prior_R_h
	f_sat ~ prior_f_sat
	
	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	

	R = LilGuys.calc_R_ell(data.xi .- d_xi, data.eta .- d_eta, 
                           ellipticity, position_angle)

	prof = LilGuys.Plummer(r_s=R_h, M=1)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R)
	L_bg_space = 1 / area_tot
	
    L = @. ((1-f_sat) * L_bg_space * data.L_bg
            + f_sat *  L_sat_space * data.L_sat)
    L = max.(L, LL_min)
    LL = sum(log.(L))
	
	Turing.@addlogprob!(LL)
end


@model function sersic_ell_model(data::GaiaData;
        area_tot::Real, R_max::Real,
        prior_tangent = prior_tangent, 
        prior_R_h = prior_R_h, 
        prior_f_sat = prior_f_sat
    )

	d_xi  ~ prior_tangent
	d_eta ~ prior_tangent
	R_h ~ prior_R_h
	f_sat ~ prior_f_sat
	n ~ Uniform(0.4, 12)

	ellipticity ~ Uniform(0, 0.99)
	position_angle ~ Uniform(0, 180)
	

	
	R = LilGuys.calc_R_ell(data.xi .- d_xi, data.eta .- d_eta,
                           ellipticity, position_angle)

	b_n = LilGuys.guess_b_n(n)
	prof = LilGuys.Sersic(n=n, R_h=R_h, _b_n=b_n)

	mtot = LilGuys.mass_2D(prof, R_max)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R) / mtot
	L_bg_space = 1 / area_tot
	
    L = @. ((1-f_sat) * L_bg_space * data.L_bg
            + f_sat *  L_sat_space * data.L_sat)
    L = max.(L, LL_min)
    LL = sum(log.(L))
	
	Turing.@addlogprob!(LL)
end


@model function sersic_model(data::GaiaData;
        area_tot::Real, R_max::Real,
        prior_tangent = prior_tangent, 
        prior_R_h = prior_R_h, 
        prior_f_sat = prior_f_sat
    )

	d_xi  ~ prior_tangent
	d_eta ~ prior_tangent
	R_h ~ prior_R_h
	f_sat ~ prior_f_sat
	n ~ Uniform(0.4, 12)


	
	R = @. sqrt((data.xi - d_xi)^2 + (data.eta - d_eta)^2)

	b_n = LilGuys.guess_b_n(n)
	prof = LilGuys.Sersic(n=n, R_h=R_h, _b_n=b_n)

	mtot = LilGuys.mass_2D(prof, R_max)
	
	L_sat_space = @. LilGuys.surface_density.(prof, R) / mtot
	L_bg_space = 1 / area_tot
	
    L = @. ((1-f_sat) * L_bg_space * data.L_bg
            + f_sat *  L_sat_space * data.L_sat)
    L = max.(L, LL_min)
    LL = sum(log.(L))
	
	
	Turing.@addlogprob!(LL)
end



"""
    summarize(samples)

Sumerizes a Turing chain. This loosely wraps Turings version
except adds median and quantile uncertainties.
"""
function summarize(samples)
	chain_summary = DataFrame(Turing.summarize(samples))
	Nvar = size(chain_summary, 1)

	medians = zeros(Nvar)
	err_lows = zeros(Nvar)
	err_highs = zeros(Nvar)
	for i in 1:Nvar
        l, m, h = quantile(samples[:, i, :], [PVALUE, 0.5, 1-PVALUE])
		medians[i] = m
		err_lows[i] = m - l
		err_highs[i] = h - m
	end

	chain_summary[!, :parameters] = string.(chain_summary.parameters)
	chain_summary[!, :median] = medians
	chain_summary[!, :lower_error] = err_lows
	chain_summary[!, :upper_error] = err_highs

	chain_summary
end


function to_measurements(A::AbstractMatrix; pvalue=PVALUE)
	m = dropdims(median(A, dims=2), dims=2)
	h = quantile.(eachrow(A), 1-pvalue)
	l = quantile.(eachrow(A), pvalue)

	return Measurement.(m, m .- l, h .- m)
end

function to_measurement(A::AbstractVector; pvalue=PVALUE)
    m = median(A)
	h = quantile(A, 1-pvalue)
	l = quantile(A, pvalue)

	return Measurement(m, m .- l, h .- m)
end

