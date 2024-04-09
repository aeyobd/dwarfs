### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# ╔═╡ 83d84568-f46b-11ee-32cb-013c7608aab0
begin
	import Pkg; Pkg.activate()
	using Plots; gr()

	import LilGuys as lguys
	import HDF5
	import CSV
	using DataFrames
	import Arya
end

# ╔═╡ b4d6d376-ad71-46ac-8767-7f3a3e0a8b0a
dirname = "/home/dboyea/sculptor/isolation/1e4/fiducial/"

# ╔═╡ 7e3528d8-3698-44bf-9a6b-562117bc36d1
cd(dirname)

# ╔═╡ c85cb861-51af-4063-ae17-9619e13cc57a
out = lguys.Output("out")

# ╔═╡ dba2b598-170d-448e-9eb6-41e30eeace1b
i = 1

# ╔═╡ d196627f-11f8-4462-8be0-136dd672eb38
begin 
	function add_centres!(filename, label=filename)
		cens = CSV.read("data/$filename", DataFrame)
		x_c = transpose(Matrix(cens[:, ["x", "y", "z"]]))
		push!(x_cs, x_c)
		push!(labels, label)
	end
	
	x_cs = []
	labels = []
	add_centres!("centres.csv", "ss")
	add_centres!("centres_com.csv", "com")
	add_centres!("centres_pot.csv", "potential")
	#add_centres!("rapha_ss_centres.csv", "rapha")
end

# ╔═╡ d68b5cc0-3e88-4dd2-8fe4-2d0991f85e06
x_cs

# ╔═╡ bf616146-46a5-427b-a776-02dc935722d3
# ╠═╡ disabled = true
#=╠═╡
begin 
	cens = CSV.read("centres.csv", DataFrame)
	x_c_1 = transpose(Matrix(cens[:, ["x", "y", "z"]]))

	cens2 = CSV.read("orbitfileCOD.dat", DataFrame, delim=raw" ")
	x_c_2 = transpose(Matrix(cens2[:, ["x", "y", "z"]]))

	cens = CSV.read("centres_phi.csv", DataFrame)
	x_c_phi = transpose(Matrix(cens[:, ["x", "y", "z"]]))
	
	cens = CSV.read("centres_com.csv", DataFrame)
	x_c_com = transpose(Matrix(cens[:, ["x", "y", "z"]]))

	cens = CSV.read("centres_ss0.3.csv", DataFrame)
	x_c_03= transpose(Matrix(cens[:, ["x", "y", "z"]]))
	
	x_cs = [x_c_1, x_c_2, x_c_phi, x_c_com, x_c_03]
	labels = ["ss", "ss_rapha", "phi", "com", "ss0.3"]

end
  ╠═╡ =#

# ╔═╡ 4d800299-3971-443d-b490-b62bef3682be
begin 
	anim = @animate for i in 1:1:length(out)
		lguys.plot_centre(out[i].positions, width=20, ms=2, alpha=0.3, color="black", label="")
		plot!(legend_position=:outerright, grid=false, 
			dpi=100,)
		for (label, x_c) in zip(labels, x_cs)
			scatter!(x_c[1, i:i], x_c[2, i:i], label=label)
		end	
	
	end
	gif(anim, "figures/isolation.gif", fps=10)
end

# ╔═╡ bf447376-68fb-4d09-846e-04a22bf4d377
begin 
	anim3 = @animate for i in 1:1:length(out)
		plot(xlims=(-20, 20), ylims=(-20, 20),
			legend_position=:outerright, grid=false, dpi=100,
		)
		scatter!(out[i].positions[1, :], out[i].positions[2, :], 
			label="", ms=2, alpha=0.3, color="black")
		plot!()
		for (label, x_c) in zip(labels, x_cs)
			scatter!(x_c[1, i:i], x_c[2, i:i], label=label)
		end	
	
	end
	gif(anim3, fps=10)
end

# ╔═╡ b134e6cb-3d63-46a9-a882-b2dc5b13a27f
# ╠═╡ disabled = true
#=╠═╡
begin 
	anim2 = @animate for i in 1:10:length(out)
		plot(legend=false, grid=false, axis=false, dpi=100, aspect_ratio=1)
		lguys.plot_centre!(out[i].positions, width=1000, ms=2, alpha=0.3)
	end
	gif(anim2, fps=1)
end
  ╠═╡ =#

# ╔═╡ e488fd62-f04f-4a67-aa5d-57a58438c13f
lguys.plot_xyz_layout(x_cs..., labels=labels, aspect_ratio=1)

# ╔═╡ Cell order:
# ╠═83d84568-f46b-11ee-32cb-013c7608aab0
# ╠═b4d6d376-ad71-46ac-8767-7f3a3e0a8b0a
# ╠═7e3528d8-3698-44bf-9a6b-562117bc36d1
# ╠═c85cb861-51af-4063-ae17-9619e13cc57a
# ╠═dba2b598-170d-448e-9eb6-41e30eeace1b
# ╠═d196627f-11f8-4462-8be0-136dd672eb38
# ╠═d68b5cc0-3e88-4dd2-8fe4-2d0991f85e06
# ╠═bf616146-46a5-427b-a776-02dc935722d3
# ╠═4d800299-3971-443d-b490-b62bef3682be
# ╠═bf447376-68fb-4d09-846e-04a22bf4d377
# ╠═b134e6cb-3d63-46a9-a882-b2dc5b13a27f
# ╠═e488fd62-f04f-4a67-aa5d-57a58438c13f
