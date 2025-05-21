using SpecialFunctions: erf
using LilGuys
using CairoMakie
import Base: -, reverse

import Agama
import Agama: py2mat, mat2py, py2vec, py2f


potential_dir = ENV["DWARFS_ROOT"] * "/agama/potentials/"

⊕(x::Real, y::Real) = sqrt(x^2 + y^2)

# vasiliev units
V_V2KMS = 1
V_M2MSUN = 232_500
V_R2KPC = 1
V_T2GYR = 0.97779

F = Float64


function load_agama_potential(filename)
    return Agama.Potential(joinpath(potential_dir, filename))
end
