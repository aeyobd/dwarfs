### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ ecfd2d72-08e0-11ef-23b5-55f91aef7d76
begin 
	using Pkg; Pkg.activate()
	using CairoMakie;
	using CSV, DataFrames

	import LilGuys as lguys

	using Arya
end

# ╔═╡ b0fa3f7c-be01-475c-b731-126226bfe337
using LinearAlgebra

# ╔═╡ f2284d88-090c-4836-a503-d230eec5569e
using Colors

# ╔═╡ 0a5520ea-604e-476e-a761-e12a0121cf27
Makie.update_theme!(px_per_unit=0.5)

# ╔═╡ 043e6913-7b92-4d55-962f-33c1d6e6b524
cd("/cosma/home/durham/dc-boye1/sculptor/orbits/orbit1")

# ╔═╡ 8ee597b8-e5b1-4d5b-a8d6-03fb1a32544e
starsfile = "../../isolation/1e6/stars/exp2d_stars.hdf5"

# ╔═╡ 44a1f33e-18ae-4a2c-9973-2db7acf666cb
let 
	using HDF5

	f = h5open(starsfile)
	p_idx = f["index"][:]
	global probabilities = f["probabilities"][:][sortperm(p_idx)]
	p_idx = sort(p_idx)
	close(f)
	
end

# ╔═╡ 88535209-6ff9-45ee-90ef-4939bd79789c
begin 
	out =  lguys.Output("out/combined.hdf5")
	
	cens = CSV.read("out/centres.csv", DataFrames.DataFrame)
	x_cen = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	v_cen = transpose(Matrix(cens[:, ["vx", "vy", "vz"]]))
	out.x_cen .= x_cen
	out.v_cen .= v_cen

	out
end

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

	hm = Arya.hist2d!(ax, snap.positions[2, :] .- x_cen[2, i], snap.positions[3, :] .- x_cen[3, i], bins = bins, weights = probabilities[snap.index], colormap=:greys,
		colorscale=log10,
		colorrange=colorrange
	)
	
	record(fig, "sculptor_stars.mp4", idxs, framerate = framerate) do i
		snap = out[i]
		x = snap.positions[2, :] .- x_cen[2, i]
		y = snap.positions[3, :] .- x_cen[3, i]
		ms = probabilities[snap.index]
		H, x_e, y_e = Arya.histogram2d(x, y, bins, weights=ms)

		println(maximum(H))
		println(minimum(H[H .> 0]))
		
		hm[3] = H
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

	hm = Arya.hist2d!(ax, snap.positions[2, :] .- x_cen[2, i], snap.positions[3, :] .- x_cen[3, i], bins = bins,  colormap=:greys,
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


# ╔═╡ 6d461b4f-7d42-4853-8288-7f8f399c39ab
Base.@kwdef struct Camera
    eye::Vector{Float64}     # Camera position
    center::Vector{Float64}  = [0,0,0]
    up::Vector{Float64}  = [0,0,1]
    fov::Float64             # Field of view in degrees
    aspect_ratio::Float64   = 1
    near_clip::Float64       = 0
    far_clip::Float64        = 1e6

    # Constructor to ensure the camera is always properly initialized
    function Camera(eye, center, up, fov, aspect_ratio, near_clip, far_clip)
        new(eye, center, up, fov, aspect_ratio, near_clip, far_clip)
    end
end

# ╔═╡ b31a51ba-c843-4ad8-9f94-6e20ccbc352f
function project(camera::Camera, point::AbstractVector{Float64})
    # Define camera basis vectors
    zaxis = normalize(camera.center - camera.eye)  # Forward
    xaxis = normalize(cross(zaxis, camera.up))  # Right
    yaxis = cross(xaxis, zaxis)  # True up

    # View matrix
    view = [xaxis yaxis -zaxis camera.eye;
            0 0 0 1]

    # Perspective projection matrix
    f = cotd(camera.fov / 2)  # Focal length
    perspective = [f/camera.aspect_ratio 0 0 0;
                   0 f 0 0;
                   0 0 (camera.far_clip + camera.near_clip) / (camera.near_clip - camera.far_clip) (2 * camera.far_clip * camera.near_clip) / (camera.near_clip - camera.far_clip);
                   0 0 -1 0]

    # Convert point to homogeneous coordinates
    point_homogeneous = [point; 1]

    # Apply view matrix
    viewed_point = view * point_homogeneous

    # Apply perspective matrix
    perspective_point = perspective * viewed_point

    # Perform perspective divide to get 2D coordinates
    if perspective_point[4] != 0
        return perspective_point[1:2] ./ perspective_point[4]
    else
        return nothing  # Return nothing if point is at infinity or invalid
    end
end


# ╔═╡ 85e8933e-49c5-4f93-82a1-e52ec5925277
function project(camera::Camera, points::AbstractMatrix)
	N = size(points, 2)
    projected_points  = Matrix{eltype(points)}(undef, 2, N )
    for i in 1:N
	    projected_points[:, i] .= project(camera, points[:, i])
    end

	return projected_points
end

# ╔═╡ a8a4196e-68cf-4385-9231-7bb21a79bf72
begin 
	struct Detector
	    xbins::AbstractVector  # x-coordinate bins
	    ybins::AbstractVector  # y-coordinate bins
	    counts::Matrix{Int}     # Matrix to store counts of projected points
	end
	# Constructor to initialize the detector with the size of the bins
function Detector(xbins::AbstractVector{T}, ybins::AbstractVector{T}) where T <: Real
	return Detector(xbins, ybins, zeros(Int, length(ybins)-1, length(xbins)-1))
end
end

# ╔═╡ a7349fee-b223-4d44-b67e-ab2dca97f840
begin 
	function register!(detector::Detector, x::Real, y::Real)
		# Find the appropriate bin (careful with bounds)
		xbin_idx = searchsortedfirst(detector.xbins,x)
		ybin_idx = searchsortedfirst(detector.ybins, y)
		if xbin_idx > 1 && xbin_idx <= length(detector.xbins)
			if ybin_idx > 1 && ybin_idx <= length(detector.ybins)
				detector.counts[ybin_idx-1, xbin_idx-1] += 1
			end
		end
	end

	function register!(detector::Detector, x::AbstractVector, y::AbstractVector)
		for (a, b) in zip(x, y)
			register!(detector, a, b)
		end
	end
end

# ╔═╡ fee4d3c0-3be9-443d-951a-6a4124e9d828
function detect_points!(camera::Camera, detector::Detector, points::Matrix{Float64})
    for point in eachcol(points)
        projected_point = project(camera, point)
        register!(detector, projected_point...)
    end
end


# ╔═╡ aebaafc7-e9ce-4546-82f8-a618a8955048
detector = Detector(LinRange(-10, 10, 300), LinRange(-10, 10, 300))

# ╔═╡ 071af567-d7cb-47b0-b5a6-5937b7ff4e79
function loglog(x) 
	if !isfinite(x) || 1>=x
		return NaN
	else
		log(log(x))
	end
end

# ╔═╡ b4e9969a-20bf-4858-9cd5-042ab9f35228
function plot_detector(detector)
	heatmap(detector.xbins, detector.ybins, detector.counts', 
		colorscale=loglog, colorrange=(2,1e4), 
		axis=(;aspect=1))
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
plot_detector(detector)

# ╔═╡ bc5d65c1-efb2-4d8f-b154-201540dd8947
LinRange(-1, 1, 1000) isa AbstractVector

# ╔═╡ 1463a550-23b3-4426-be93-5bbbd2b0c98f


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

# ╔═╡ 7777bd75-cddf-4a3d-8bc3-27f59e38b1db
cam = Camera(eye=[1000,0,20], fov=25)

# ╔═╡ 8769428a-c954-4a92-8786-58d5923d3b47
project(cam, [0.,0.,0.])

# ╔═╡ 3aeccc9b-0608-4217-af4e-dcc41f12d2e1
detect_points!(cam, detector, out[end].positions)

# ╔═╡ e179fb4f-330e-4ade-9341-3e17da14d05d
project(cam, snap_i.positions)

# ╔═╡ 40c18d1e-2245-4834-af48-3a9999d90cb2
lguys.get_y(snap_i.positions)

# ╔═╡ 6d7ab3e9-7f64-4ad9-8342-c3c3252d8037
p_idx = 1:1000:length(snap_i)

# ╔═╡ d0910d76-3226-4774-8b87-de5cc97c1a49
typeof(p_idx)

# ╔═╡ e8d3f8a3-ffb7-40d0-b5ad-1c9d54bf3d25
supertypes(StepRange)

# ╔═╡ a98977b9-96a0-4f81-b601-fda3bba8731a
function Makie.Point3f(snap::lguys.Snapshot)
	return Point3f[r for r in eachcol(snap_i.positions)]
end

# ╔═╡ f196ab04-d343-45c9-9e61-ff167ae7c5ba
positions = Point3f(snap_i[p_idx])

# ╔═╡ 3e59070e-f84f-4beb-97e4-d172d4a8087d
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ 7e6e74ec-278d-4c52-989b-81457ed210f4
# ╠═╡ disabled = true
#=╠═╡
begin
	scene = Scene()

	α_inv = 1

	p = scatter!(scene, positions;
		markersize=3, alpha=1. /α_inv, color="white"
	)

	scene.backgroundcolor = colorant"black"

	cam3d!(scene, fov=20)

	move_camera!(scene, 0)

	scene
end
  ╠═╡ =#

# ╔═╡ 32fdbee4-b292-4a8d-aef5-57a6ebe35742
function move_camera!(scene, time; R_0 = 1000, z_0=50, ϕ0=π/3, ω=2π)
	ϕ = ϕ0 + ω * time
	eyeposition = (R_0*cos(ϕ), R_0*sin(ϕ), z_0)
	lookat = (0,0,0)
	update_cam!(scene, eyeposition, lookat)
end

# ╔═╡ b5e359bc-3388-4869-b868-1d7f8d8c8f5c
# ╠═╡ disabled = true
#=╠═╡
let
	
	scene = Scene()

	α_inv = 100

	snap = out[100]

	x = lguys.get_x(snap)
	println(size(x))
	y = lguys.get_y(snap)
	z = lguys.get_z(snap)

	filt = 1:10:length(snap)

	positions = Observable(snap.positions)
	p = scatter!(scene, positions;
		markersize=3, alpha=1/α_inv, color="white"
	)


	scene.backgroundcolor = colorant"black"

	cam3d!(scene, fov=20)

	R_0 = 1000
	z_0 = 50
	ϕ0 = π/3
		
	eyeposition = (R_0*cos(ϕ0), R_0*sin(ϕ0), z_0)
	lookat = (0,0,0)
	update_cam!(scene, eyeposition, lookat)

	framerate = 3
	idxs = vcat(1:10:91, 92:100:length(out))

	ω = 2π

		
	record(scene, "sculptor.mp4", eachindex(idxs), framerate = framerate) do i
		snap = out[idxs[i]]
		pos = Observable[Point{3, Float32}(x[1], x[2], x[3]) for x in eachcol(snap.positions)]

		#p[1] = pos
		

		ϕ = ϕ0 + ω * (i / length(idxs))
		eyeposition = (R_0*cos(ϕ), R_0*sin(ϕ), z_0)

		update_cam!(scene, eyeposition, lookat)
	end
	scene
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═ecfd2d72-08e0-11ef-23b5-55f91aef7d76
# ╠═0a5520ea-604e-476e-a761-e12a0121cf27
# ╠═043e6913-7b92-4d55-962f-33c1d6e6b524
# ╠═8ee597b8-e5b1-4d5b-a8d6-03fb1a32544e
# ╠═44a1f33e-18ae-4a2c-9973-2db7acf666cb
# ╠═b0fa3f7c-be01-475c-b731-126226bfe337
# ╠═88535209-6ff9-45ee-90ef-4939bd79789c
# ╠═57497573-5d20-4d26-8b5a-7af72580c59b
# ╠═a3858097-f307-438e-84ff-6c7be376cb8e
# ╠═c05c2e3c-192f-479a-9b4c-9a0558275a5d
# ╠═e449d00f-dda2-4bf7-81aa-fe4454310ec9
# ╟─18d8144c-f1f6-494e-90a0-527c72f3880b
# ╠═a5914e73-5866-48c8-a98e-ea0c3fb53178
# ╠═6d461b4f-7d42-4853-8288-7f8f399c39ab
# ╠═b31a51ba-c843-4ad8-9f94-6e20ccbc352f
# ╠═8769428a-c954-4a92-8786-58d5923d3b47
# ╠═85e8933e-49c5-4f93-82a1-e52ec5925277
# ╠═a8a4196e-68cf-4385-9231-7bb21a79bf72
# ╠═a7349fee-b223-4d44-b67e-ab2dca97f840
# ╠═fee4d3c0-3be9-443d-951a-6a4124e9d828
# ╠═b4e9969a-20bf-4858-9cd5-042ab9f35228
# ╠═aebaafc7-e9ce-4546-82f8-a618a8955048
# ╠═3aeccc9b-0608-4217-af4e-dcc41f12d2e1
# ╠═071af567-d7cb-47b0-b5a6-5937b7ff4e79
# ╠═3664612b-4106-4b3c-bed1-a17d6f4a7a34
# ╠═0d252adf-f414-451e-b09c-ab6c79ae30dd
# ╠═b99c72c7-d049-4ed8-9c4c-535fe863d176
# ╠═966618dc-85c0-41ef-9306-7998904ffc5d
# ╠═56a5585b-305d-4a9b-ae56-a0bfc7d6117e
# ╠═4aa1cf3d-cb47-4b8f-90d9-a9aa88296946
# ╠═bc5d65c1-efb2-4d8f-b154-201540dd8947
# ╠═1463a550-23b3-4426-be93-5bbbd2b0c98f
# ╠═751c780d-0d52-446e-a53f-fcf535819f6c
# ╠═7777bd75-cddf-4a3d-8bc3-27f59e38b1db
# ╠═e179fb4f-330e-4ade-9341-3e17da14d05d
# ╠═40c18d1e-2245-4834-af48-3a9999d90cb2
# ╠═f2284d88-090c-4836-a503-d230eec5569e
# ╠═6d7ab3e9-7f64-4ad9-8342-c3c3252d8037
# ╠═d0910d76-3226-4774-8b87-de5cc97c1a49
# ╠═e8d3f8a3-ffb7-40d0-b5ad-1c9d54bf3d25
# ╠═f196ab04-d343-45c9-9e61-ff167ae7c5ba
# ╠═a98977b9-96a0-4f81-b601-fda3bba8731a
# ╠═3e59070e-f84f-4beb-97e4-d172d4a8087d
# ╠═7e6e74ec-278d-4c52-989b-81457ed210f4
# ╠═32fdbee4-b292-4a8d-aef5-57a6ebe35742
# ╠═b5e359bc-3388-4869-b868-1d7f8d8c8f5c
