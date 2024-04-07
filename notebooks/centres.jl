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
dirname = "/home/dboyea/sculptor/isolation/1e4"

# ╔═╡ 7e3528d8-3698-44bf-9a6b-562117bc36d1
cd(dirname)

# ╔═╡ c85cb861-51af-4063-ae17-9619e13cc57a
out = lguys.Output("out")

# ╔═╡ dba2b598-170d-448e-9eb6-41e30eeace1b
i = 1

# ╔═╡ bf616146-46a5-427b-a776-02dc935722d3
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

# ╔═╡ 44c654c5-d70b-4895-b079-50c89c19918c
# ╠═╡ skip_as_script = true
#=╠═╡
begin 
	lguys.plot_centre(out[i].positions, width=10, ms=2, alpha=0.2)
	scatter!(x_c_1[1, i:i], x_c_1[2, i:i], label="ss")
	scatter!(x_c_2[1, i:i], x_c_2[2, i:i], label="ss rapha")

end
  ╠═╡ =#

# ╔═╡ 4d800299-3971-443d-b490-b62bef3682be
begin 
	anim = @animate for i in 1:1:length(out)
		lguys.plot_centre(out[i].positions, width=10, ms=2, alpha=0.3, color="black", label="")
		plot!(legend_position=:outerright, grid=false, 
			dpi=100,)
		for (label, x_c) in zip(labels, x_cs)
			scatter!(x_c[1, i:i], x_c[2, i:i], label=label)
		end	
	
	end
	gif(anim, "isolation.gif", fps=10)
end

# ╔═╡ b134e6cb-3d63-46a9-a882-b2dc5b13a27f
begin 
	anim2 = @animate for i in 1:10:length(out)
		plot(legend=false, grid=false, axis=false, dpi=100, aspect_ratio=1)
		lguys.plot_centre!(out[i].positions, width=1000, ms=2, alpha=0.3)
	end
	gif(anim2, fps=1)
end

# ╔═╡ 82c9502c-7678-4f4f-a7c6-81f0859a5b6c
cens2

# ╔═╡ 5f48a4b3-d19d-4dcd-9323-d7d6df35bc22
lguys.plot_xyz(x_cs..., labels=labels)[1]

# ╔═╡ 68690dae-6695-4cd9-b109-d58b266a09a5
lguys.plot_xyz(x_cs..., labels=labels)[2]

# ╔═╡ 782498c9-ddc9-4190-9163-a7a41326b08e
lguys.plot_xyz(x_cs..., labels=labels)[3]

# ╔═╡ Cell order:
# ╠═83d84568-f46b-11ee-32cb-013c7608aab0
# ╠═b4d6d376-ad71-46ac-8767-7f3a3e0a8b0a
# ╠═7e3528d8-3698-44bf-9a6b-562117bc36d1
# ╠═c85cb861-51af-4063-ae17-9619e13cc57a
# ╠═dba2b598-170d-448e-9eb6-41e30eeace1b
# ╠═44c654c5-d70b-4895-b079-50c89c19918c
# ╠═bf616146-46a5-427b-a776-02dc935722d3
# ╠═4d800299-3971-443d-b490-b62bef3682be
# ╠═b134e6cb-3d63-46a9-a882-b2dc5b13a27f
# ╠═82c9502c-7678-4f4f-a7c6-81f0859a5b6c
# ╠═5f48a4b3-d19d-4dcd-9323-d7d6df35bc22
# ╠═68690dae-6695-4cd9-b109-d58b266a09a5
# ╠═782498c9-ddc9-4190-9163-a7a41326b08e
