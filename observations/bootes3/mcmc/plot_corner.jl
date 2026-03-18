using Arya
using CSV
using DataFrames
using PairPlots
using LilGuys
using CairoMakie
CairoMakie.activate!(type=:png)

FIGDIR = "./figures"



function (@main)(ARGS)
    @assert length(ARGS) == 1
    inname = ARGS[1]
    samples = CSV.read(inname, DataFrame)
    summary = CSV.read(replace(inname, "samples"=>"summary"), DataFrame)

    cols = summary.parameters

    for col in ["R_h", "R_s", "R_s_outer"]
        if col in cols
            samples[!, "log_$col"] = log10.(samples[!, col])
            cols = cols[cols .!= [col]]
            push!(cols, "log_$col")
        end
    end


    name, _ = replace(splitext(inname), "samples."=>"")

    @savefig "corner.$name" pairplot(samples[:, cols])

    return 0
end
