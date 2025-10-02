using CSV, DataFrames
using Printf
import StatsBase: median, quantile

suffix = ARGS[1]
df_both = CSV.read(suffix, DataFrame)

θs_both = mod1.(df_both.Θ_grad, 360.) .- 360

θ_m_both= median(θs_both)
θ_both_err = quantile(θs_both, [0.16, 0.84]) .- median(θs_both)

r_grad_m_both = median(df_both.r_grad)
r_grad_both_err = quantile(df_both.r_grad, [0.16, 0.84]) .- r_grad_m_both

dsigma = median(df_both.dlσ_dlR)
dsigma_err = quantile(df_both.dlσ_dlR, [0.16, 0.84]) .- dsigma

@printf "theta = %0.1f + %0.1f - %0.1f " θ_m_both θ_both_err[2] θ_both_err[1]
println()
@printf "r = %0.2f + %0.2f - %0.2f " r_grad_m_both r_grad_both_err[2] r_grad_both_err[1] 
println()
@printf "dsigma = %0.2f + %0.2f - %0.2f " dsigma dsigma_err[2] dsigma_err[1] 
println()
