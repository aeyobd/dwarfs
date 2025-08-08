using TOML
using LilGuys
using Printf


function (@main)(ARGS)
    filename = ARGS[1]

    df = TOML.parsefile(filename)


    for key in ["ra", "dec", "distance", "pmra", "pmdec", "radial_velocity"]
        @printf "%s = %s\n" key df[key]
    end

    xi = df["x_i"], df["y_i"], df["z_i"]

    println()


    @printf "t_i = %0.2f\n" df["t_i"] * T2GYR
    @printf "pericentre = %0.2f\n" df["pericentre"]
    @printf "apocentre = %0.2f\n" df["apocentre"]
    @printf "t last peri = %0.2f\n" df["time_last_peri"] * T2GYR


    @printf "x_i = [%0.2f %0.2f %0.2f]\n" xi...


    vi = df["v_x_i"], df["v_y_i"], df["v_z_i"]
    vi = vi .* V2KMS

    @printf "v_i = [%0.2f %0.2f %0.2f]\n" vi...

end
