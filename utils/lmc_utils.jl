using PythonCall
agama = pyimport("agama")
np = pyimport("numpy")


# some simple python utilities
py2f(x) = pyconvert(Float64, x)
py2vec(x) = pyconvert(Vector{Float64}, x)
py2mat(x) = pyconvert(Matrix{Float64}, x)'



"""
    make_lmc_pot(Mlmc, r_s; kwargs...)

Creates a LMC agama potential as truncated NFW similar to vasiliev24. 
Mlmc is the mass (in our code units) and r_s is the scale radius (also in our code units).
kwargs passed to agama.Potential.
"""
function make_lmc_pot(Mlmc, r_s; kwargs...)
    Φ_lmc = agama.Potential(;
            type="Spheroid", alpha=1, beta=3, gamma=1, 
            scaleradius = r_s,
            mass = Mlmc, 
            outercutoffradius = 10r_s, 
            cutoffStrength=4,
            kwargs...
        )
    return Φ_lmc
end


"""
    calc_lmc_orbit(pot, gc_coord, kwargs...)

Computes an LMC orbit in the given potential with the present-day position `gc_coord`.
Returns an `Orbit` object.

# Arguments
- `Mlmc::function` a function which returns the mass of the LMC (our code units) given the time
- `reflex_motion::Bool=true` whether to include reflex motion in the calculation
- `dynamical_friction::Bool=true` whether to include (Chandrasakar dynamical friction) in the calculation
- `Λ = nothing`. The value of the argument of the Coloumb logarithm, 
controlling the scaling of Dynamical friction
- `σv`. A function returning the local halo velocity dispersion given distance from origin
- `vasiliev_units::Bool = false` If true, assumes the potential `pot` is in the alternate unit system in vasiliev++
- `time` Time to integrate to. If negative, integrates backwards in time. Code units.
- `r_s`
- `timestep`
Additional kwargs passed to `leap_frog`

"""
function calc_lmc_orbit(pot, gc_coord; 
        Mlmc = t->15,
        reflex_motion = true,
        dynamical_friction = true,
        Λ = nothing,
        σv = nothing,
        units = :vasiliev,
        time = -10 / T2GYR, 
        r_s = 8.5 * (Mlmc / 10)^0.6,
        timestep = :adaptive,
        kwargs...
    )

    
    if units == :vasiliev
        acc_scale = (V_V2KMS/ V2KMS)^2  * (V_R2KPC / R2KPC)^-1
        v_scale = V2KMS / V_V2KMS
        m_scale = M2MSUN / V_M2MSUN
    elseif units == :code
        v_scale = m_scale = acc_scale = 1
    end
    

    if σv === nothing
        calc_σv_interp(pot, log_r=LinRange(2, -2.0, 1000))
    end

    ρ(x) = py2f(pot.density(x))
    
    Φ_lmc = make_lmc_pot(Mlmc(0) * m_scale, r_s)


    f_fric(pos, vel, t) = dynamical_friction * a_dyn_friction(pos, vel, 
        r_s=r_s,
        σv=σv, ρ=ρ, M=Mlmc(t) * m_scale, Λ = sqrt(calc_r(pos) / 0.8r_s )
    )


    f_grav(pos, vel, t) = py2vec(
        pot.force(pos) - reflex_motion * Mlmc(t) / Mlmc(0) * Φ_lmc.force(-pos)
    )

    f_acc(pos, vel, t) = acc_scale * (f_grav(pos, vel*v_scale, t) .+ f_fric(pos, vel *v_scale, t) )

    orbit = leap_frog(gc_coord, f_acc; time=time, timestep=timestep, kwargs...)

    return orbit
end


"""
    make_lmc_mw_pot_from_orbit

Given the MW potential and the orbit of the LMC,
returns an agama.Potential object representing the combined potential.
"""
function make_lmc_mw_pot_from_orbit(pot, orbit;
        reflex_motion=true,
        Mlmc=t->15,
        r_s=8.5 , 
        units = :vasiliev,     
    )
    
    time = orbit.time

    if units == :vasiliev
        m_scale = M2MSUN / V_M2MSUN
        time *= T2GYR / V_T2GYR
    elseif units == :code
        m_scale = 1
    end
    
    if orbit.time[2] < orbit.time[1]
        position = orbit.position[:, end:-1:1]
        time = reverse(time)
    end

    centre = vcat(time', position)

    centre = PyArray(centre')
    scale = PyArray(vcat(time',  Mlmc.(time)' ./ Mlmc(0), ones(length(time))')')

    pot_lmc = make_lmc_pot(m_scale * Mlmc(0), r_s, center=centre,scale=scale )

    
    if reflex_motion != 0
        N = length(time)
        a_reflex = pot_lmc.force(zeros(length(time), 3), t=time)
        nptime = PyMatrix(reshape(time,(:, 1)))
        mat_reflex = np.hstack([nptime, -a_reflex])        
        pot_reflex = agama.Potential(type="UniformAcceleration", file=mat_reflex)
        
        Φ = agama.Potential(pot, pot_lmc, pot_reflex)
    else
        Φ = agama.Potential(pot, pot_lmc)
    end

    return Φ
end



"""
    make_lmc_mw_pot(pot, lmc_gc; kwargs...)

Given the MW potential and the LMC initial conditions, computes the orbit 
of the LMC and return the combined potential and the orbit.
Arguments are as in `calc_lmc_orbit`
"""
function make_lmc_mw_pot(pot, lmc_gc;
        reflex_motion=true,
        Mlmc=t->15,
        r_s=8.5 * (Mlmc / 10)^0.6, 
        units = :vasiliev, 
        kwargs...
    
    )
    orbit = calc_lmc_orbit(pot, lmc_gc; Mlmc=Mlmc, r_s=r_s, units=units, reflex_motion=reflex_motion, kwargs...)
    Φ = make_lmc_mw_pot_from_orbit(pot, orbit; Mlmc=Mlmc, r_s=r_s, units=units, reflex_motion=reflex_motion)

    return Φ, orbit
end



function calc_σv_func(pot_halo)
    gridr = np.logspace(1, 3, 16)

    df = agama.DistributionFunction(type="quasispherical", potential=pot_halo)
    gm = agama.GalaxyModel(pot_halo, df)

    x = np.column_stack((gridr, gridr*0, gridr*0))
    
    sigmas = py2mat(gm.moments(x, dens=false, vel=false))[1, :] .^ 0.5

    println("sigmas = ", sigmas)
    sigmafnc = agama.Spline(gridr, sigmas)
    return sigmafnc
end


"""
    calc_σv_interp(pot; log_r = LinRange(-3, 5, 100000))

Computes the velocity dispersion as a function of radius by integrating the Jeans equation.

"""
function calc_σv_interp(pot; log_r = LinRange(-3, 5, 100000))
    x0 = [1/√2, 0., 1/√2] # direction

    if !issorted(log_r)
        @warn "log_r is not sorted. Sorting."
        log_r = sort(log_r)
    end

    log_r = reverse(log_r)
    radii = 10 .^ log_r

    positions = x0' .* radii
    acc = pyconvert(Matrix{Float64}, pot.force(positions))
    rho = pyconvert(Array{Float64}, pot.density(positions))

    a = calc_r(acc')
    
    dr = abs.(LilGuys.gradient(radii))

    σ2 = 1 ./ rho .* cumsum(rho .* a .* dr)

    return LilGuys.lerp(radii, sqrt.(σ2))
end
