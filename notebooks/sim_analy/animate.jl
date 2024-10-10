### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ ecfd2d72-08e0-11ef-23b5-55f91aef7d76
begin 
	using Pkg; Pkg.activate()
	using GLMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ 222d659a-fe83-4069-ac9b-80ca52a56dd5
using Metal

# ╔═╡ b0fa3f7c-be01-475c-b731-126226bfe337
using LinearAlgebra

# ╔═╡ 0a5520ea-604e-476e-a761-e12a0121cf27
Makie.update_theme!(px_per_unit=0.5)

# ╔═╡ 8ee597b8-e5b1-4d5b-a8d6-03fb1a32544e
starsfile = "../../isolation/1e6/stars/exp2d_stars.hdf5"

# ╔═╡ ce3a5c49-0160-474b-bd9e-6cca0d5079a1
model_dir = "/Users/daniel/dwarfs/analysis/sculptor/1e4_V31_r3.2/orbit_mean/"

# ╔═╡ 88535209-6ff9-45ee-90ef-4939bd79789c
out =  lguys.Output(model_dir)

# ╔═╡ a3858097-f307-438e-84ff-6c7be376cb8e
let
	fig = Figure()
	r_max = 10
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc",
	limits=(-r_max, r_max, -r_max, r_max))

	bins = LinRange(-r_max, r_max, 1000)
	colorrange=(1e-11, 10e-3)

	framerate = 10
	idxs = vcat(1:1:91, 92:10:length(out))

	i = 1
	snap = out[i]

	hm = Arya.hist2d!(ax, snap.positions[2, :] .- snap.x_cen[2], snap.positions[3, :] .- snap.x_cen[3], bins = bins, weights = snap.weights[snap.index], colormap=:greys,
		colorscale=log10,
		colorrange=colorrange
	)
	
	record(fig, "sculptor_stars.mp4", idxs, framerate = framerate) do i
		snap = out[i]
		x = snap.positions[2, :] .- snap.x_cen[2]
		y = snap.positions[3, :] .- snap.x_cen[3]
		#ms = snap[snap.index]
		H = Arya.histogram2d(x, y, bins, weights=snap.weights)

		println(maximum(H))
		println(minimum(H[H .> 0]))
		
		hm[3] = H.values
	end

	fig
end

# ╔═╡ c05c2e3c-192f-479a-9b4c-9a0558275a5d
let
	fig = Figure()
	r_max = 10
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc",
	limits=(-r_max, r_max, -r_max, r_max))

	bins = LinRange(-r_max, r_max, 1000)
	colorrange=(1e-11, 10e-3)

	framerate = 10
	idxs = vcat(1:1:91, 92:10:length(out))

	i = 1
	snap = out[i]

	hm = Arya.hist2d!(ax, snap.positions[2, :] .- snap.x_cen[2], snap.positions[3, :] .- snap.x_cen[3], bins = bins,  colormap=:greys,
		colorscale=log10,
		colorrange=colorrange
	)
	
	# record(fig, "sculptor_stars.mp4", idxs, framerate = framerate) do i
	# 	snap = out[i]
	# 	x = snap.positions[2, :] .- x_cen[2, i]
	# 	y = snap.positions[3, :] .- x_cen[3, i]
	# 	H, x_e, y_e = Arya.histogram2d(x, y, bins)

	# 	println(maximum(H))
	# 	println(minimum(H[H .> 0]))
		
	# 	hm[3] = H
	# end

	fig
end

# ╔═╡ e449d00f-dda2-4bf7-81aa-fe4454310ec9


# ╔═╡ 18d8144c-f1f6-494e-90a0-527c72f3880b
md"""
# 3d animation
much too slow without a GPU  :(
"""

# ╔═╡ a5914e73-5866-48c8-a98e-ea0c3fb53178
# Define the Camera struct
Base.@kwdef struct Camera
    eye::Vector{Float64}     # Camera position
    center::Vector{Float64}  = [0, 0, 0]
    up::Vector{Float64}      = [0, 0, 1]
    fov::Float64             # Field of view in degrees
    aspect_ratio::Float64    = 1
    near_clip::Float64       = 0.1
    far_clip::Float64        = 1e6
end

# ╔═╡ b04bf8ef-d948-4472-9c77-e3df72184f50
"""
    compute_view_matrix(camera::Camera) -> Matrix{Float64}

Computes the view matrix for the given `camera`. The view matrix transforms 
points from world coordinates to camera coordinates.

# Arguments:
- `camera::Camera`: The camera struct defining the position, target, and up direction of the camera.

# Returns:
- A 4x4 view matrix for transforming points into camera space.
"""
function compute_view_matrix(camera::Camera)
    z_axis = normalize(camera.eye - camera.center)  # Forward direction
    x_axis = normalize(cross(camera.up, z_axis))    # Right direction
    y_axis = cross(z_axis, x_axis)                 # True up direction

    # View matrix: transformation from world space to camera space
    view = [
        x_axis[1] y_axis[1] z_axis[1] 0.0;
        x_axis[2] y_axis[2] z_axis[2] 0.0;
        x_axis[3] y_axis[3] z_axis[3] 0.0;
        -dot(x_axis, camera.eye) -dot(y_axis, camera.eye) -dot(z_axis, camera.eye) 1.0
    ]
    return view'
end

# ╔═╡ f7c1f411-89a9-4ba0-bbe8-7c730835240d
"""
    compute_projection_matrix(camera::Camera) -> Matrix{Float64}

Computes the perspective projection matrix for the given `camera`. This matrix 
projects 3D points into 2D normalized device coordinates (NDC).

# Arguments:
- `camera::Camera`: The camera struct defining the field of view, aspect ratio, 
  and near/far clipping planes.

# Returns:
- A 4x4 projection matrix for projecting points into 2D.
"""
function compute_projection_matrix(camera::Camera)
    fov_rad = deg2rad(camera.fov)
    f = 1.0 / tan(fov_rad / 2.0)  # Focal length based on fov

    # Perspective projection matrix
    proj = [
        f / camera.aspect_ratio 0 0 0;
        0 f 0 0;
        0 0 (camera.far_clip + camera.near_clip) / (camera.near_clip - camera.far_clip) -1;
        0 0 2.0 * camera.far_clip * camera.near_clip / (camera.near_clip - camera.far_clip) 0
    ]
    return proj
end

# ╔═╡ 6d461b4f-7d42-4853-8288-7f8f399c39ab
"""
    project_points(camera::Camera, points::Matrix{Float64}) -> Matrix{Float64}

Projects 3D points into 2D screen space using the given `camera` parameters. 
First, it applies the view transformation to transform points into the camera's 
coordinate system. Then, it applies the projection matrix to project points 
into normalized device coordinates (NDC), and finally scales them into 2D screen space.

# Arguments:
- `camera::Camera`: The camera struct that defines the position, orientation, 
  and projection properties.
- `points::Matrix{Float64}`: A matrix of 3D points of size `N x 3` where each 
  row represents a point `[x, y, z]`.

# Returns:
- A matrix of projected 2D points of size `N x 2`, where each row represents 
  the projected point `[x', y']`.
"""
function project_points(camera::Camera, points::AbstractMatrix{<:Real})
    # Compute view and projection matrices
    view_matrix = compute_view_matrix(camera)
    proj_matrix = compute_projection_matrix(camera)

	N = size(points, 2)

	points_4d = vcat(points, ones(eltype(points), N)')

	if points isa MtlArray{<:Real}
		view_matrix = MtlArray{Float32}(view_matrix)
		proj_matrix = MtlArray{Float32}(proj_matrix)
	end
	
    # Transform the points into camera space (view transformation)
    camera_space_points = view_matrix * points_4d

    # Apply perspective projection (project into 2D screen space)
    projected_points = proj_matrix * camera_space_points
    #projected_points = camera_space_points

    # Perform perspective divide to normalize the points
    ndc_points = projected_points[1:3, :] ./ projected_points[4:4, :]

    # Scale to screen coordinates ([-1, 1] range in NDC to screen)
    screen_points = ndc_points[1:2, :]'  
    return screen_points
end

# ╔═╡ 071af567-d7cb-47b0-b5a6-5937b7ff4e79
function loglog(x) 
	if !isfinite(x) || 1>=x
		return NaN
	else
		log(log(x))
	end
end

# ╔═╡ 3664612b-4106-4b3c-bed1-a17d6f4a7a34
expexp(x) = exp(exp(x))

# ╔═╡ 0d252adf-f414-451e-b09c-ab6c79ae30dd
Makie.inverse_transform(::typeof(loglog)) = expexp

# ╔═╡ b99c72c7-d049-4ed8-9c4c-535fe863d176
Makie.defaultlimits(::typeof(loglog)) = (1., 5.)

# ╔═╡ 966618dc-85c0-41ef-9306-7998904ffc5d
Makie.defined_interval(::typeof(loglog)) = Makie.IntervalSets.OpenInterval(1, Inf)

# ╔═╡ 56a5585b-305d-4a9b-ae56-a0bfc7d6117e
loglog(0.1)

# ╔═╡ 4aa1cf3d-cb47-4b8f-90d9-a9aa88296946
let 

	cam = Camera(eye= 300 .* [0, 1,0], fov=5)
	xy = project_points(cam, out[end].positions)

	fig, ax = FigAxis(
		limits = 100 .* (-1, 1, -1, 1),
		xgridvisible=false, 
		ygridvisible=false,

	)
	println(xy[1:10, :])
	scatter!(xy[:, 1], xy[:, 2], alpha=0.5, markersize=3)

	fig
end

# ╔═╡ 751c780d-0d52-446e-a53f-fcf535819f6c
snap_i = out[1]

# ╔═╡ 57497573-5d20-4d26-8b5a-7af72580c59b
let
	fig = Figure()
	r_max = 150
	ax = Axis(fig[1,1], aspect=1,
	xlabel = "y / kpc", ylabel="z / kpc", title="dark matter",
	limits=(-r_max, r_max, -r_max, r_max))

	bins = LinRange(-r_max, r_max, 300)
	colorrange=(1, 1e3)

	framerate = 10
	idxs = vcat(1:91, 92:10:length(out))

	hm = Arya.hist2d!(ax, snap_i.positions[2, :], snap_i.positions[3, :], bins = bins, colorscale=log10, colorrange=colorrange, colormap=:greys)
	
	# record(fig, "sculptor.mp4", idxs, framerate = framerate) do i
	# 	snap = out[i]
	# 	x = snap.positions[2, :]
	# 	y = snap.positions[3, :]
	# 	H, x_e, y_e = Arya.histogram2d(x, y, bins)

	# 	println(maximum(H))
	# 	println(minimum(H[H .> 0]))
		
	# 	hm[3] = H
	# end

	fig
end

# ╔═╡ 3aeccc9b-0608-4217-af4e-dcc41f12d2e1
vcat(snap_i.positions, ones(10000)')

# ╔═╡ 39a61f89-bb98-4425-90bd-d01175680cf9
vcat(snap_i.positions, ones(10000)')

# ╔═╡ 40c18d1e-2245-4834-af48-3a9999d90cb2
lguys.get_y(snap_i.positions)

# ╔═╡ a98977b9-96a0-4f81-b601-fda3bba8731a
function Makie.Point3f(snap::lguys.Snapshot)
	return Point3f[r for r in eachcol(snap_i.positions)]
end

# ╔═╡ 32fdbee4-b292-4a8d-aef5-57a6ebe35742
function move_camera!(scene, time; R_0 = 1000, z_0=50, ϕ0=π/3, ω=2π)
	ϕ = ϕ0 + ω * time
	eyeposition = (R_0*cos(ϕ), R_0*sin(ϕ), z_0)
	lookat = (0,0,0)
	update_cam!(scene, eyeposition, lookat)
end

# ╔═╡ 7e6e74ec-278d-4c52-989b-81457ed210f4
begin
	scene = Scene()

	α_inv = 1

	p = scatter!(scene, snap_i.positions;
		markersize=3, alpha=1. /α_inv, color="white"
	)

	scene.backgroundcolor = colorant"black"

	cam3d!(scene, fov=20)

	move_camera!(scene, 0)

	scene
end

# ╔═╡ 3e59070e-f84f-4beb-97e4-d172d4a8087d
function make_scene(positions)
	scene = Scene()
	scene.backgroundcolor = colorant"black"

	cam3d!(scene, fov=20)

	move_camera!(scene, 0)

	p = scatter!(scene, positions;
		markersize=3, alpha=1. /α_inv, color="white"
	)

	return scene
end

# ╔═╡ b5e359bc-3388-4869-b868-1d7f8d8c8f5c
# ╠═╡ disabled = true
#=╠═╡
let

	idxs = 1:1:length(out)

	
	scene = Scene()

	α_inv = 10

	snap = out[idxs[1]]

	x = lguys.get_x(snap)
	println(size(x))
	y = lguys.get_y(snap)
	z = lguys.get_z(snap)

	positions = Observable(snap.positions)
	p = scatter!(scene, positions;
		markersize=3, alpha=1/α_inv, color="white"
	)


	scene.backgroundcolor = colorant"black"

	cam3d!(scene, fov=20)

	move_camera!(scene, 0)


	framerate = 60

	ω = 0.5
	

		
	record(scene, "sculptor.mp4", idxs, framerate = framerate) do i
		snap = out[i]
		
		positions[] = snap.positions
		t = out.times[i] * lguys.T2GYR

		
		move_camera!(scene, t, ω=ω)
	end
	scene
end
  ╠═╡ =#

# ╔═╡ c5ce0e70-0530-4789-885e-81b98d3d7378
let

	idxs = 1:100:length(out)

	fov = 5
	
	fig = Figure()

	ax = Axis(fig[1, 1],
		limits=100 .*(-1, 1, -1, 1),
		aspect=DataAspect(),
	)

	hidedecorations!(ax)

	
	α_inv = 10
	ω = 0.0

	
	snap = out[idxs[1]]

	r0 = 300
	ϕ = π/2

	t0 = π/2
	eye(t) = [sin(ϕ) * sin(ω*t - t0), sin(ϕ) * cos(ω*t - t0), cos(ϕ), ]
	
	cam = Camera(eye= r0 .* eye(0), fov=fov)

	xy = project_points(cam, out[1].positions)
	positions = Observable(xy)
	
	p = scatter!(positions, color=:white, alpha=0.05, markersize=4)
	
	ax.backgroundcolor = colorant"black"

	framerate = 60

	

		
	record(fig, "sculptor.mp4", idxs, framerate = framerate) do i
		snap = out[i]
		t = out.times[i] * lguys.T2GYR


		cam = Camera(eye= r0 .* eye(t), fov=fov)
		
		positions[] = project_points(cam, snap.positions)

	end
	fig
end

# ╔═╡ 4f59f9d2-00d3-47c4-86d1-8d08ef6373db
cam = Camera(eye=[300, 0, 0], fov=0.5)

# ╔═╡ 0733ede2-4275-4502-98a3-fdead289d2a8
@time for i in eachindex(out)
	xy = project_points(cam, out[i].positions)
end

# ╔═╡ 92552458-ccf6-45ab-a1b4-1ecb5a7b27b7
@time for i in eachindex(out)
	xy = project_points(cam, out[i].positions)
	Arya.histogram2d(xy[1, :], xy[2, :], 1000, limits=100 .*(-2, 2, -2, 2))
end

# ╔═╡ f372bae2-2402-4fe7-b72d-b73a7805b280
@time for i in eachindex(out)
	xy = project_points(cam, MtlArray{Float32}(out[i].positions))
end

# ╔═╡ Cell order:
# ╠═ecfd2d72-08e0-11ef-23b5-55f91aef7d76
# ╠═222d659a-fe83-4069-ac9b-80ca52a56dd5
# ╠═0a5520ea-604e-476e-a761-e12a0121cf27
# ╠═8ee597b8-e5b1-4d5b-a8d6-03fb1a32544e
# ╠═b0fa3f7c-be01-475c-b731-126226bfe337
# ╠═ce3a5c49-0160-474b-bd9e-6cca0d5079a1
# ╠═88535209-6ff9-45ee-90ef-4939bd79789c
# ╠═57497573-5d20-4d26-8b5a-7af72580c59b
# ╠═a3858097-f307-438e-84ff-6c7be376cb8e
# ╠═c05c2e3c-192f-479a-9b4c-9a0558275a5d
# ╠═e449d00f-dda2-4bf7-81aa-fe4454310ec9
# ╟─18d8144c-f1f6-494e-90a0-527c72f3880b
# ╠═a5914e73-5866-48c8-a98e-ea0c3fb53178
# ╠═b04bf8ef-d948-4472-9c77-e3df72184f50
# ╠═f7c1f411-89a9-4ba0-bbe8-7c730835240d
# ╠═6d461b4f-7d42-4853-8288-7f8f399c39ab
# ╠═3aeccc9b-0608-4217-af4e-dcc41f12d2e1
# ╠═071af567-d7cb-47b0-b5a6-5937b7ff4e79
# ╠═3664612b-4106-4b3c-bed1-a17d6f4a7a34
# ╠═0d252adf-f414-451e-b09c-ab6c79ae30dd
# ╠═b99c72c7-d049-4ed8-9c4c-535fe863d176
# ╠═966618dc-85c0-41ef-9306-7998904ffc5d
# ╠═56a5585b-305d-4a9b-ae56-a0bfc7d6117e
# ╠═4aa1cf3d-cb47-4b8f-90d9-a9aa88296946
# ╠═39a61f89-bb98-4425-90bd-d01175680cf9
# ╠═751c780d-0d52-446e-a53f-fcf535819f6c
# ╠═40c18d1e-2245-4834-af48-3a9999d90cb2
# ╠═a98977b9-96a0-4f81-b601-fda3bba8731a
# ╠═3e59070e-f84f-4beb-97e4-d172d4a8087d
# ╠═7e6e74ec-278d-4c52-989b-81457ed210f4
# ╠═32fdbee4-b292-4a8d-aef5-57a6ebe35742
# ╠═b5e359bc-3388-4869-b868-1d7f8d8c8f5c
# ╠═c5ce0e70-0530-4789-885e-81b98d3d7378
# ╠═4f59f9d2-00d3-47c4-86d1-8d08ef6373db
# ╠═0733ede2-4275-4502-98a3-fdead289d2a8
# ╠═92552458-ccf6-45ab-a1b4-1ecb5a7b27b7
# ╠═f372bae2-2402-4fe7-b72d-b73a7805b280
