function ellipse!(radius, ellipticity, position_angle; kwargs...)
	t = LinRange(0, 2π, 1000) 

	h = 1/sqrt(1 - ellipticity)
	a = radius * h
	b = radius / h

	x = @. a * cos(t) 
	y = @. b * sin(t)

	θ = deg2rad(position_angle)
	x1 = @. x * sin(θ) - y*cos(θ)
	y1 = @. x * cos(θ) + y*sin(θ)

	rs = LilGuys.calc_R_ell(x1, y1, ellipticity, position_angle)

	@assert all(rs .≈ radius)
	lines!(x1, y1; kwargs...)
end
