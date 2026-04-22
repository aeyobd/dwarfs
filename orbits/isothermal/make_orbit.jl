using Agama
using LilGuys
pot = Agama.Potential(file="agama_potential.ini")
Φ(r) = Agama.potential(pot, [0, 0, r])

function peri_apo_to_E_L(peri, apo)
	L = sqrt(
		(Φ(peri) - Φ(apo)) / 
			(1/2 * (1/apo^2 - 1/peri^2))
	)

	E = Φ(peri) + L^2 / 2peri^2

	@assert Φ(apo) + L^2 / 2apo^2 ≈ E
	return E, L
end



function orbital_period(peri, apo)
	E, L = peri_apo_to_E_L(peri, apo)
	return 2*LilGuys.integrate(
		r -> 1 / sqrt(2*(E - Φ(r)) - L^2 / r^2),
		peri * (1 + 1e-12), apo*(1 - 1e-12)
	)
end



function (@main)(ARGS)
    T_orb = orbital_period(peri, apo)
    orbit = LilGuys.agama_orbit(pot, Galactocentric([0, 0, apo], [0, L/apo * V2KMS, 0]), timerange=(0, 10*T_orb))
    
end

