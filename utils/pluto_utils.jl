using PlutoUI
using CairoMakie

SIMDIR = joinpath(ENV["DWARFS_SIMS"])
DWARFSDIR = joinpath(ENV["DWARFS_ROOT"])

function notebook_inputs(; kwargs...)
	return PlutoUI.combine() do Child
		
		user_inputs = [
			md""" $(string(name)): $(
				Child(name, obj)
			)"""
			
			for (name, obj) in kwargs
		]
		
		md"""
		#### Inputs
		$(user_inputs)
		"""
	end
end


CairoMakie.activate!(type=:png)
