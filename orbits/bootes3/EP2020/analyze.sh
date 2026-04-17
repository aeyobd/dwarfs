julia ../../derive_quantiles.jl --galaxy bootes3  --p-value 0.022750131948179205 . # 2 sigma
julia ../../derive_quantiles.jl --galaxy bootes3 . derived_quantiles_3sigma.toml
julia ../../derive_quantiles.jl --galaxy bootes3 --p-value 0.15865525393145696 . derived_quantiles_1sigma.toml 
