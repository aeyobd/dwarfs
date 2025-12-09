using LilGuys

probs_inner = LilGuys.read_hdf5_table("../exp2d_rs0.10/probabilities_stars.hdf5")
probs_outer = LilGuys.read_hdf5_table("../exp2d_rs0.21/probabilities_stars.hdf5")

f_outer = 0.315

@assert probs_inner.index == probs_outer.index
@assert sum(probs_inner.probability) ≈ 1
@assert sum(probs_outer.probability) ≈ 1

probs_inner.probability = probs_inner.probability * (1-f_outer) .+ f_outer * probs_outer.probability


@assert sum(probs_inner.probability) ≈ 1

LilGuys.write_hdf5_table("probabilities_stars.hdf5", probs_inner)


