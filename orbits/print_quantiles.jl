
import TOML


function (@main)(ARGS)
    allprops = TOML.parsefile(ARGS[1])

    for (orbit, props) in allprops
        println("[$orbit]")
        for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
            val = props[key * "_median"]
            println(key, " = ", round(val, digits=2))
        end
        println
    end
end
