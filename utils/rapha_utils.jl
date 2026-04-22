using LilGuys


const R2r = LilGuys.r_h(LilGuys.Exp2D()) / LilGuys.R_h(LilGuys.Exp2D())
const α_exp = LilGuys.R_h(LilGuys.Exp2D())
const α_exp_cusp = LilGuys.r_circ_max(LilGuys.ExpCusp())
const α_nfw = LilGuys.r_circ_max(LilGuys.NFW(r_s=1, M_s=1))

integrate = LilGuys.integrate
lerp = LilGuys.lerp
find_zero = LilGuys.find_zero


⊕(x, y) = sqrt(x^2 + y^2)

Base.@kwdef struct TidalNFW <: LilGuys.GeneralNFW
    r_circ_max::Real
    v_circ_max::Real
    r_cut::Real
    kappa::Real = 0.3
    _nfw::LilGuys.NFW = LilGuys.NFW(r_circ_max=r_circ_max, v_circ_max=v_circ_max)
    r_s::Real = r_circ_max / α_nfw
end


function LilGuys.density(h::TidalNFW, r)
    return LilGuys.density(h._nfw, r) * exp(-r / h.r_cut) / (1 + h._nfw.r_s / h.r_cut)^h.kappa
end


function r_v_rel(M_max_rel)
    r_rel = find_zero(r -> LilGuys.v_circ_EN21(r)^2 * r - M_max_rel, 0.5)
    v_rel = LilGuys.v_circ_EN21(r_rel)
    return r_rel, v_rel
end


function final_halo(r_max, v_max, M_max_rel)
    # Eq. 9 in EN2021
    r_cut = r_max * 0.44 * M_max_rel^0.44 * (1 - M_max_rel^0.3)^(-1.1)
    r_rel, v_rel = r_v_rel(M_max_rel)

    h =  TidalNFW(r_circ_max=r_max, v_circ_max=v_max, r_cut=r_cut)

    r_scale = r_max*r_rel / LilGuys.r_circ_max(h)
    v_scale = v_max * v_rel / LilGuys.v_circ_max(h) 

    return TidalNFW(r_circ_max=r_max*r_scale, v_circ_max=v_max*v_scale, r_cut=r_cut * r_scale)
end


function calc_σv(halo, R_h=R_h)
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



function density_from_f(f, Ψ, r)
	Ψ_0 = Ψ(r)
	# Eq 4.43 in BT2008
	
	return 4π * LilGuys.integrate(
		E -> f(E) * sqrt(2*(Ψ_0 - E)),
		0, Ψ_0
	)
end

function make_stellar_density(f_star, Ψ, r)
	ρ =  density_from_f.([f_star], Ψ, r)

	ρ[end-1:end] .= 0
	# ρ[1] = density_from_df.(df_stars, Φ, 0)
	return InterpProfile(LilGuys.lerp(r, ρ))
end



function density_of_states(Ψ, E)
	r_max = LilGuys.find_zero(r -> Ψ(r) + -E, (1e-5, 1e8))
	

	@assert 1 >= E > 0
	(4π)^2 * LilGuys.integrate(r-> r^2 * sqrt(2 * (-E + Ψ(r))), 
					   0, r_max)
end


struct InterpProfile <: LilGuys.SphericalProfile
	interp
end

LilGuys.density(pot::InterpProfile, r) = pot.interp(r)


function calc_r_h(ρ; R_end=Inf)
	
	M(r) =  4π * LilGuys.integrate(r -> ρ(r) * r^2, 0, r)

	Mtot = M(R_end)
	return LilGuys.find_zero(r -> M(r) - 1/2 * Mtot, (1e-8, R_end))
end



function calc_mass(ρ, r)
	4π * LilGuys.integrate(r -> ρ(r) * r^2, 0, r)
end



function final_ρ_star(initial_props, M_max_rel)
    @info "calculating dfs for M_rel = $M_max_rel"
    halo_f = final_halo(initial_props.r_max_i, initial_props.v_max_i, M_max_rel)
	Ψ_f(r) = -LilGuys.potential(halo_f, r)
    ϵ_max_f = Ψ_f(initial_props.r_i[1] * 0.5) 

    @info "final energies"
	x_f, dNde_star_f = final_stellar_energies(initial_props.x_i, initial_props.dNde_star_i, M_max_rel)
	ϵ_f = (1 .- x_f) .* (ϵ_max_f)
	g_f = density_of_states.(Ψ_f, ϵ_f)

	f_star_f = dNde_star_f ./ g_f
    r_f = [find_zero(r -> Ψ_f(r) - ϵ, (1e-8, 1e8)) for (ϵ, r_i) in zip(ϵ_f, initial_props.r_i)]

   return make_stellar_density(LilGuys.lerp(ϵ_f, f_star_f), Ψ_f, r_f)
end


function final_ρ_star(halo_i::LilGuys.SphericalProfile, M_max_rel; 
        kwargs...
    )

    initial_props = find_initial_df(halo_i; kwargs...)
    return final_ρ_star(initial_props, M_max_rel)
end


function find_initial_df(halo_i; 
        r_range = (1e-4, 1e4), N_r = 200,
        α_star = 3, β_star=3, E_s
    )
    r_dm = logrange(r_range[1], r_range[2], N_r) |> collect
	Ψ_i(r) = -LilGuys.potential(halo_i, r)
    ϵ_i = Ψ_i.(r_dm)

    ϵ_max_i = Ψ_i(r_dm[1] * 0.5)
    x_i = 1 .- ϵ_i / ϵ_max_i

    dNde_star_i = dN_dE_poly.(x_i, E_s, α_star, β_star)

    r_max_i = LilGuys.r_circ_max(halo_i)
    v_max_i = LilGuys.v_circ_max(halo_i)

    return (;
        r_max_i = r_max_i,
        v_max_i = v_max_i,
        dNde_star_i = dNde_star_i,
        x_i = x_i,
        r_i = r_dm,
    )
end


function find_initial_df(halo_i, prof_stars; r_range = (1e-4, 1e4), N_r = 100,

    )
    r_dm = logrange(r_range[1], r_range[2], N_r) |> collect
	Ψ_i(r) = -LilGuys.potential(halo_i, r)
    ϵ_i = Ψ_i.(r_dm)

    ϵ_max_i = Ψ_i(r_dm[1] * 0.5)
    x_i = 1 .- ϵ_i / ϵ_max_i

    g_i = density_of_states.(Ψ_i, ϵ_i)

    f_star_i = LilGuys.DistributionFunction(
        LilGuys.density.(prof_stars, r_dm), 
        Ψ_i.(r_dm), r_dm
    ).(ϵ_i)
	dNde_star_i = f_star_i .* g_i


    r_max_i = LilGuys.r_circ_max(halo_i)
    v_max_i = LilGuys.v_circ_max(halo_i)

    return (;
        r_max_i = r_max_i,
        v_max_i = v_max_i,
        dNde_star_i = dNde_star_i,
        x_i = x_i,
        r_i = r_dm,
    )
end



function calc_r_h(halo_i, E_s, α_star=3, β_star = 3; r_range=(1e-4, 1e4), N_r=100)
	r_dm = logrange(r_range[1], r_range[2], N_r) |> collect
	Ψ_i(r) = -LilGuys.potential(halo_i, r)
    ϵ_i = Ψ_i.(r_dm)

    ϵ_max_i = Ψ_i(r_range[1] * 0.5) # 
    x_i = 1 .- ϵ_i / ϵ_max_i


    dNde_star_i = dN_dE_poly.(x_i, E_s, α_star, β_star)
    g_i = density_of_states.(Ψ_i, ϵ_i)
	
	f_star_i = dNde_star_i ./ g_i

	ρ = make_stellar_density(LilGuys.lerp(ϵ_i, f_star_i), Ψ_i, r_dm)
	return calc_r_h(ρ.interp, R_end=1e4)
end

