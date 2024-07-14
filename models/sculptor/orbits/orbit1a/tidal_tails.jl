### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ fb8bb8ba-34ad-11ef-23e6-1d890b60e0b9
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ 20ac4d7c-d835-4a7f-9600-261e43f4b290
using FITSIO

# ╔═╡ 9c7035e7-c1e7-40d5-8ab6-38f0bb682111
md"""
# Tidal Tails
A detailed analysis of the stars in sculptor
"""

# ╔═╡ cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
starsfile = "stars/exp2d_rs0.13_today.fits"

# ╔═╡ 7a92c896-7552-4f35-9761-5709d23e9adf
FITS(starsfile, "r") do f
	global stars = DataFrame(f[2])
end

# ╔═╡ 8dbc941f-287e-49f4-8bcf-123e3851f015
begin 
	cens = CSV.read(joinpath(".", "out/centres.csv"), DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
end

# ╔═╡ 8a08fc3b-f112-46df-b4d7-10c22fc2f8eb
snap_cen = lguys.Snapshot(x_cen, v_cen, ones(size(x_cen, 2)))

# ╔═╡ a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
obs_today_filename = "../../mc_orbits/orbit1.toml"

# ╔═╡ c4008b83-b61b-4baa-9fd5-9cced2dc6be8
import TOML

# ╔═╡ e37559b2-229c-4a37-b516-c6cb7c022b71
orbit_props = TOML.parsefile(joinpath(".", "orbital_properties.toml"))

# ╔═╡ 0ba05fdb-f859-4381-b2d0-145aa04f7bbf
obs_today_file = TOML.parsefile(obs_today_filename)

# ╔═╡ 4ef955fb-a813-46ad-8f71-7c8a3d371eee
obs_today_icrs = lguys.ICRS(;
	ra=obs_today_file["ra"], dec=obs_today_file["dec"],
	distance=obs_today_file["distance"],
	pm_ra=obs_today_file["pm_ra"],
	pm_dec=obs_today_file["pm_dec"],
	radial_velocity=obs_today_file["radial_velocity"],
)

# ╔═╡ c3a3129e-18c6-4348-81c0-b03c5836b785
frame = lguys.HelioRest

# ╔═╡ 075ae901-bbbd-4d10-9d91-c393fc86a8e7
function make_sample(snap; 
	cen=nothing, rel_p_cut=1e-15, r_max=Inf,
	Frame=frame
)
	snap_stars = snap

	
	obs_pred = lguys.to_sky(snap_stars, SkyFrame=Frame)
	
	obs_df = DataFrame(; 
	collect(key => [getproperty(o, key) for o in obs_pred]
			for key in [:ra, :dec, :pm_ra, :pm_dec, :distance, :radial_velocity])...
		
	)

	obs_df[!, "index"] = snap_stars.index
	obs_df[!, "probability"] = ones(length(snap_stars))

	if cen !== nothing
		obs_c_galcen = lguys.Galactocentric(x=cen.x, y=cen.y, z=cen.z, 
			v_x=cen.v_x, v_y=cen.v_y, v_z=cen.v_z)
		obs_c = lguys.transform(lguys.ICRS, obs_c_galcen)
	
		cen_df = DataFrame(ra=obs_c.ra, dec=obs_c.dec, pm_ra=obs_c.pm_ra, pm_dec=obs_c.pm_dec, radial_velocity=obs_c.radial_velocity, distance=obs_c.distance, index=-1, probability=0.0)
	
		obs_df = append!(cen_df, obs_df)
		
		obs_df[!, "xi"], obs_df[!, "eta"] = lguys.to_tangent(obs_df.ra, obs_df.dec, obs_c.ra, obs_c.dec)
		obs_df[!, "r_ell"] = @. 60 * sqrt(obs_df.xi^2 + obs_df.eta^2)

	end
	
	rename!(obs_df, "pm_ra"=>"pmra")
	rename!(obs_df, "pm_dec"=>"pmdec")

	return obs_df
end

# ╔═╡ 61d3487b-2fad-4b8b-b46b-ddc0c42f61e9
sky_orbit = make_sample(snap_cen, Frame=frame)

# ╔═╡ 910267ee-aa39-4f04-961b-70f0862d27e2
obs_today = lguys.transform(frame, obs_today_icrs)

# ╔═╡ 4c1d6f6f-e257-4126-9b4a-8e5aa8470295
function ra_dec_axis(ddeg=5; kwargs...)
	fig = Figure(;kwargs...)
	
	dy = ddeg
	dx = dy * 1/cosd(obs_today.dec)
	limits = (obs_today.ra .+ (-dx, dx), obs_today.dec .+ (-dy, dy))

	ax = Axis(fig[1,1],
		xlabel="RA / degrees",
		ylabel="dec / degrees",
		limits=limits,
		aspect = 1,
		xgridvisible=false,
		ygridvisible=false
		
	)

	return fig, ax
end

# ╔═╡ 6d0cff99-3efb-4405-a11d-f13200fa5334
idx_f = orbit_props["idx_f"]

# ╔═╡ 816b9db9-26c6-4ac8-9a46-82209d2cdc85
idx_orbit = idx_f - 180: idx_f + 80


# ╔═╡ e9e35643-168e-4e87-a880-6831b46145c7
let 
	fig, ax = ra_dec_axis()

	bins = 100
	limits = ax.limits.val
	x = stars.ra
	y = stars.dec


	hi = Arya.histogram2d(x, y, bins, weights=stars.probability, limits=limits)
	areas = diff(hi.xbins) .* (diff(hi.ybins)')
	hi.values ./= areas
	
		
	h = heatmap!(hi, colorscale=log10, colorrange=(1e-10, maximum(hi.values)))

	lines!(sky_orbit.ra[idx_orbit], sky_orbit.dec[idx_orbit])
	
	Colorbar(fig[1, 2], h,
		label="stellar density"
	)

	fig
end

# ╔═╡ e15f2d0d-3d4f-4772-8f3e-a00f2360cddc
function spherical_to_cartesian(ra::Vector{Float64}, dec::Vector{Float64})
    x = cosd.(dec) .* cosd.(ra)
    y = cosd.(dec) .* sind.(ra)
    z = sind.(dec)
    return hcat(x, y, z)'
end

# ╔═╡ 70627dd3-7b6f-43fb-8ec5-7a4ab7cfe1e8
import LinearAlgebra: svd, norm, cross, dot, I, inv

# ╔═╡ 2c194883-a76e-48c1-8cea-76e0ae860045
# Function to fit a great circle to given Cartesian coordinates
function fit_great_circle(coords::AbstractMatrix{Float64})
    # Perform Singular Value Decomposition
    U, S, V = svd(coords')
    # The normal vector to the best-fit plane is the right singular vector corresponding to the smallest singular value
    normal_vector = V[:, 3]
    return normal_vector
end


# ╔═╡ 7e36537e-0d1a-409e-bd2a-985ff6b05955
function fit_great_circle(ra, dec)
	coords = spherical_to_cartesian(ra, dec)
	return fit_great_circle(coords)
end

# ╔═╡ ee17890b-d84b-4212-ab4f-e9da73c503d6
# Function to rotate coordinates using a rotation matrix
function rotate_coordinates(coords::AbstractMatrix{Float64}, rotation_matrix::AbstractMatrix{Float64})
    return rotation_matrix * coords
end

# ╔═╡ a7dfc2cf-db9f-4c6b-b3a8-85ee8810445c
# Function to create a rotation matrix that aligns normal_vector with the z-axis
function create_rotation_matrix(normal_vector::Vector{Float64})
    # Normalize the normal vector
    normal_vector /= norm(normal_vector)
    # Calculate the axis and angle of rotation
    z_axis = [0.0, 0.0, 1.0]
    axis = cross(normal_vector, z_axis)
    if norm(axis) == 0
        return I(3)  # Identity matrix, no rotation needed
    end
    axis /= norm(axis)
    angle = acos(dot(normal_vector, z_axis))
    # Create the rotation matrix using the Rodrigues' rotation formula
    K = [0.0 -axis[3] axis[2]; axis[3] 0.0 -axis[1]; -axis[2] axis[1] 0.0]
    R = I(3) + sin(angle) * K + (1 - cos(angle)) * (K * K)
    return R
end


# ╔═╡ 93556a3e-f1c9-40cb-bccf-12442464b82c
# Function to create a rotation matrix about the z-axis
function create_z_rotation_matrix(angle::Float64)
    return [cos(angle) -sin(angle) 0.0;
            sin(angle)  cos(angle) 0.0;
            0.0        0.0        1.0]
end

# ╔═╡ bad0a95e-5f8f-4556-98e3-10cea99ef93b
# Function to convert Cartesian coordinates to RA and Dec
function cartesian_to_spherical(coords::Matrix{Float64})
    x, y, z = coords[1, :], coords[2, :], coords[3, :]
    ra = atand.(y, x)
    dec = atand.(z, sqrt.(x.^2 + y.^2))
    return ra, dec
end

# ╔═╡ aeb5faa4-8eb7-49e7-8434-f1e391af99c6
coords = spherical_to_cartesian(stars.ra, stars.dec)

# ╔═╡ c0153932-83ac-45df-b9f2-2c74939741b1
lguys.unit_vector(stars.ra, stars.dec)

# ╔═╡ 3abaf91f-72cc-4fcd-8e40-c5141b83a33c
lguys.cartesian_to_sky(coords)

# ╔═╡ ba69849f-03f3-4950-8849-3ee9ab3b1a15
coords_o = spherical_to_cartesian(sky_orbit.ra[idx_orbit], sky_orbit.dec[idx_orbit])

# ╔═╡ a875cabe-0dd8-440d-8a39-dc95792d43a7
norm_vec = fit_great_circle(coords_o)

# ╔═╡ a4106b33-ab30-4065-ae53-107937a2a9a1
rp, dp = cartesian_to_spherical(rotate_coordinates(coords, create_rotation_matrix(norm_vec)))

# ╔═╡ 77a095ac-e51f-4628-bc3c-917966bf334e
norm_vec

# ╔═╡ 20cee4d0-5237-4699-93fd-77af324cfa23
sum(coords_o .* norm_vec, dims=1)

# ╔═╡ a892cee1-58b6-429f-990f-c9909d85a3a9
rpo, dpo = cartesian_to_spherical(
	rotate_coordinates(coords_o
		, create_rotation_matrix(norm_vec)))

# ╔═╡ 90036e32-cb65-4a3c-89b7-7bba3bdc836a
let
	fig, ax = FigAxis()
	scatter!(rp, dp)

	lines!(rpo, dpo, color=:red)
	fig
end

# ╔═╡ 0d80ceac-7fdf-4db1-aa36-fd8187c7df5c
ts = LinRange(0, 2π, 1000)

# ╔═╡ a180ee1c-6e94-41d5-b7db-a51b5f529e8b
circ_eq = [cos.(ts) sin.(ts) 0ts]'

# ╔═╡ eeb855c1-94f8-4b81-ba65-143bf5cfc1f2
mat = create_rotation_matrix(norm_vec)

# ╔═╡ ce046918-dd6b-4172-9517-8b465ea7a65b
mat * norm_vec

# ╔═╡ 4959e6fe-2dc5-46f3-83cd-3d0b2ba72350
mat_i = inv(mat)

# ╔═╡ 481f17e0-7a84-4b28-a6b3-6495eead1f41
lguys.plot_xyz(coords_o, mat_i * circ_eq)

# ╔═╡ Cell order:
# ╠═9c7035e7-c1e7-40d5-8ab6-38f0bb682111
# ╠═fb8bb8ba-34ad-11ef-23e6-1d890b60e0b9
# ╠═cfe54fc2-0c12-44cd-a6be-5f6cae93f68d
# ╠═20ac4d7c-d835-4a7f-9600-261e43f4b290
# ╠═7a92c896-7552-4f35-9761-5709d23e9adf
# ╠═8dbc941f-287e-49f4-8bcf-123e3851f015
# ╠═8a08fc3b-f112-46df-b4d7-10c22fc2f8eb
# ╠═075ae901-bbbd-4d10-9d91-c393fc86a8e7
# ╠═61d3487b-2fad-4b8b-b46b-ddc0c42f61e9
# ╠═a1b48fb9-af21-49e0-ae78-7a1e51c50bc4
# ╠═c4008b83-b61b-4baa-9fd5-9cced2dc6be8
# ╠═e37559b2-229c-4a37-b516-c6cb7c022b71
# ╠═0ba05fdb-f859-4381-b2d0-145aa04f7bbf
# ╠═4ef955fb-a813-46ad-8f71-7c8a3d371eee
# ╠═c3a3129e-18c6-4348-81c0-b03c5836b785
# ╠═910267ee-aa39-4f04-961b-70f0862d27e2
# ╠═4c1d6f6f-e257-4126-9b4a-8e5aa8470295
# ╠═6d0cff99-3efb-4405-a11d-f13200fa5334
# ╠═e9e35643-168e-4e87-a880-6831b46145c7
# ╠═816b9db9-26c6-4ac8-9a46-82209d2cdc85
# ╠═e15f2d0d-3d4f-4772-8f3e-a00f2360cddc
# ╠═70627dd3-7b6f-43fb-8ec5-7a4ab7cfe1e8
# ╠═2c194883-a76e-48c1-8cea-76e0ae860045
# ╠═7e36537e-0d1a-409e-bd2a-985ff6b05955
# ╠═ee17890b-d84b-4212-ab4f-e9da73c503d6
# ╠═a7dfc2cf-db9f-4c6b-b3a8-85ee8810445c
# ╠═93556a3e-f1c9-40cb-bccf-12442464b82c
# ╠═a875cabe-0dd8-440d-8a39-dc95792d43a7
# ╠═bad0a95e-5f8f-4556-98e3-10cea99ef93b
# ╠═aeb5faa4-8eb7-49e7-8434-f1e391af99c6
# ╠═c0153932-83ac-45df-b9f2-2c74939741b1
# ╠═3abaf91f-72cc-4fcd-8e40-c5141b83a33c
# ╠═20cee4d0-5237-4699-93fd-77af324cfa23
# ╠═a4106b33-ab30-4065-ae53-107937a2a9a1
# ╠═77a095ac-e51f-4628-bc3c-917966bf334e
# ╠═ba69849f-03f3-4950-8849-3ee9ab3b1a15
# ╠═a892cee1-58b6-429f-990f-c9909d85a3a9
# ╠═90036e32-cb65-4a3c-89b7-7bba3bdc836a
# ╠═0d80ceac-7fdf-4db1-aa36-fd8187c7df5c
# ╠═a180ee1c-6e94-41d5-b7db-a51b5f529e8b
# ╠═eeb855c1-94f8-4b81-ba65-143bf5cfc1f2
# ╠═ce046918-dd6b-4172-9517-8b465ea7a65b
# ╠═4959e6fe-2dc5-46f3-83cd-3d0b2ba72350
# ╠═481f17e0-7a84-4b28-a6b3-6495eead1f41
