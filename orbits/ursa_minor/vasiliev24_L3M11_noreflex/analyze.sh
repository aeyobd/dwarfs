#julia derive_quantiles.jl --galaxy ursa_minor --key pericentre_lmc . median_properties_lmc_bound.toml --p-value 0.016

julia near_quantile.jl --galaxy ursa_minor --key pericentre_lmc . selected_orbits.toml --p-value 0.016
exit 0

julia ../../derive_quantiles.jl --galaxy ursa_minor . 
julia ../../derive_quantiles.jl --galaxy ursa_minor --key pericentre_lmc . median_properties_lmc.toml 
julia ../../derive_quantiles.jl --galaxy ursa_minor --key pericentre_lmc . median_properties_lmc_2sigma.toml --p-value 0.02275013194817921
