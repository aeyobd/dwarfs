using LilGuys


const R2r = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())
const α_exp = LilGuys.R_h(LilGuys.Exp2D())
const α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())
integrate = LilGuys.integrate
lerp = LilGuys.lerp


⊕(x, y) = sqrt(x^2 + y^2)

function ExpCusp(r_max, v_max)
	r_s = r_max / α_exp_cusp
	v0 = LilGuys.v_circ_max(LilGuys.ExpCusp(1, r_s))
	M = @. (v_max/v0)^2
	halo = LilGuys.ExpCusp(M, r_s)
end


function calc_σv(r_max, v_max, R_h=R_h)
	halo = ExpCusp(r_max, v_max)
	prof_stars = LilGuys.Exp2D(R_s=R_h / α_exp)
	
	LilGuys.σv_star_mean(halo, prof_stars)
end




function σv_wolf(halo, R_h)
	return LilGuys.v_circ(halo, R_h*R2r)/sqrt(3)
end


function t_max_to_r_max_rel(t_max_rel) 
    # no 2π since both are relative quantities...
    f(r_max_rel) = r_max_rel / LilGuys.v_circ_EN21(r_max_rel) - t_max_rel

    LilGuys.find_zero(f, (1e-10, 1))
end


function ecc_factor(peri, apo)
	f_ecc = (2*apo/peri / (apo/peri + 1))^3.2
end

function rapha_final_halo(r_max, v_max, peri, apo, n = 0; V0=1, n_scale=1)
    f_ecc = ecc_factor(peri, apo)
	t_ecc = n * n_scale / f_ecc

	T_max_0 = 2π*r_max / v_max
	T_peri = 2π * peri / V0

	if T_max_0 / T_peri > 2/3 # heavy mass loss
        
        # discussion for Eq. 10
		T_asy = 0.22 * T_peri
		τ_asy = 0.65 # units of T_orb, but divided by n anyways
		
        # Eq. 12
		Y_0 = (T_max_0 - T_asy) / T_peri
        # Eq. 13
		τ = τ_asy / Y_0

        # Eq. 14
		η = 1 - exp(-2.5Y_0)
	else
        # Eq. 16
		T_asy = T_max_0  / (1 + T_max_0/T_peri)^2.2
        # Eq. 17
		τ = 1.2 * (T_max_0/T_peri)^-0.5 
        
        # discussion for Eq. 15
		η = 0.67 
		Y_0 = (T_max_0 - T_asy) / T_peri
	
	end

    # Eq. 12, 15
	Y = Y_0 * ( 1 + (t_ecc/τ)^η )^(-1/η)
    # Inverting Eq. 10
	T_max = T_peri * Y + T_asy

	r_max_rel = t_max_to_r_max_rel(T_max / T_max_0)
	v_max_rel = LilGuys.v_circ_EN21(r_max_rel)

	return r_max_rel * r_max, v_max_rel * v_max
end




function normalize_dN_dE(E, dN_dE)
    dN_dE_interp = lerp(E, dN_dE)
    return integrate(dN_dE_interp, extrema(E)...)
end

function dNde_map(x, x_max_t)
    # Eq. 9 in Errani+22
    return @. 1 / (1 + (0.85 * x/x_max_t)^12)
end

function x_t_map(x_i, x_max_t)
    # Eq. 12 in Errani+22
    x_f = @. ( 1 + (0.8 * x_i / x_max_t)^(-3) )^(1/-3)
end

function log10normal(x, μ, σ)
    z = log10(x / μ) / σ
    1 / (x * σ * √(2π) * log(10)) * exp(-z^2/2)
end

function final_stellar_energies(x_i, dN_dE_i, M_max_rel)
    # equation 10
    x_max_t = 0.77 * M_max_rel^0.43 

    dN_dE_t = @. dN_dE_i * dNde_map(x_i, x_max_t)
    x_f = x_t_map.(x_i, x_max_t)

    dN_dE_t_interp = lerp(x_i, dN_dE_t)

    # average over dNdE within 0.03dex of x
    dN_dE_conv(x_f_p) = integrate(
        x_i_pp -> dN_dE_t_interp(x_i_pp) * log10normal(x_f_p, x_t_map(x_i_pp, x_max_t), 0.03),
        minimum(x_i), 1
    )

    dN_dE_f = dN_dE_conv.(x_f)

    x_f, dN_dE_f
end


function dN_dE_poly(E::Real, E_s::Real, α::Real, β::Real)
    if 0 <= E < 1
        E^α * exp(-(E/E_s)^β)
    else
        0.0
    end
end

