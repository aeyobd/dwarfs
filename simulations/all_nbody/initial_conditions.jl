using CSV, DataFrames
using LilGuys
using PyFITS


function get_obs_props()
    df = CSV.read(joinpath(ENV["DWARFS_ROOT"], "observations/observed_properties_complete.csv"), DataFrame)

    df[!, :distance_modulus] = LilGuys.kpc2dm.(df.distance)
    df[!, :ra_err] .= 0
    df[!, :dec_err] .= 0

    return df[df.galaxyname .!= ["lmc"], :]
    df
end


function sample_initial_conditions(obs_props)
	pos_i = Point[]
	vel_i = Point[]
	halos = []

    σ_fattahi = 0.037
    σ_ludlow = 0.1
    σ_M_L = 0.2
	
	for row in eachrow(obs_props)
		icrs = LilGuys.rand_coord(Dict((k=>row[k] for k in names(row))))
		gc_i = LilGuys.transform(Galactocentric, icrs)

        M_L_star = 2.0 * 10^(σ_M_L * rand())
		L = LilGuys.mag_to_L(row.Mv)
		Mstar = L / M2MSUN * M_L_star

        vcircmax = LilGuys.vel_from_M_s_fattahi(Mstar * 10^(σ_fattahi*randn()))

        if row.galaxyname == "smc"
            # di teodoro et al,...
            vcircmax = (56 + randn() * 5) / V2KMS
        elseif row.galaxyname == "sagittarius"
            # vasiliev+2021
            vcircmax = (13.5 + randn() * 1) / V2KMS
        end

		rcircmax = LilGuys.Ludlow.solve_rmax(vcircmax, σ_ludlow * randn())
		halo = LilGuys.NFW(v_circ_max=vcircmax, r_circ_max=rcircmax)

		push!(pos_i, Point(LilGuys.position(gc_i)))
		push!(vel_i, Point(LilGuys.velocity(gc_i) / V2KMS))
		push!(halos, halo)
	end

	return pos_i, vel_i, halos
end


function read_initial_conditions(filename, halo_type=NFW; halo_attrs=[:r_s, :M_s])
    df = read_fits(filename)
    pos = [NBody.Point(r.x, r.y, r.z) for r in eachrow(df)]
    vel = [NBody.Point(r.v_x, r.v_y, r.v_z) for r in eachrow(df)]
    halos = [halo_type(; (attr=>r[attr] for attr in halo_attrs)...) for r in eachrow(df)]


    return pos, vel, halos
end


function write_initial_conditions(filename, ics)
    pos, vel, halos = ics

    df = LilGuys.to_frame(ics)
    df[!, :x] = [p[1] for p in pos]
    df[!, :y] = [p[2] for p in pos]
    df[!, :z] = [p[3] for p in pos]
    df[!, :v_x] = [p[1] for p in vel]
    df[!, :v_y] = [p[2] for p in vel]
    df[!, :v_z] = [p[3] for p in vel]

    write_fits(filename, df)
end
