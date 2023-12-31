using PyCall
using Printf


const u = PyNULL()
const astropy_coords = PyNULL()

const galcen_frame = PyNULL()
const geocen_frame = PyNULL()

function __init__()
    copy!(u, pyimport("astropy.units"))
    copy!(astropy_coords, pyimport("astropy.coordinates"))
    copy!(galcen_frame, astropy_coords.Galactocentric(
        galcen_distance = 8.29*u.kpc,
        z_sun=0*u.pc,
        galcen_v_sun = [11.1, 240.3+12.24, 7.25] * (u.km/u.s)
       ))
    copy!(geocen_frame, astropy_coords.ICRS())
end


Base.@kwdef struct Observation
    ra::F
    dec::F
    distance::F
    pm_ra::F
    pm_dec::F
    radial_velocity::F
end

Base.@kwdef struct FuzzyObservation
    ra::F
    dec::F
    distance::F
    pm_ra::F
    pm_dec::F
    radial_velocity::F

    distance_err::F
    pm_ra_err::F
    pm_dec_err::F
    radial_velocity_err::F
end





function to_galcen(obs::Observation)
    sc = astropy_coords.SkyCoord(
        ra = obs.ra * u.degree,
        dec = obs.dec * u.degree,
        distance = obs.distance * u.kpc,
        radial_velocity = obs.radial_velocity * u.km/u.s,
        pm_ra_cosdec = obs.pm_ra * u.mas/u.yr,
        pm_dec = obs.pm_dec * u.mas/u.yr
        )  
    tc = sc.transform_to(galcen_frame)

    pos = [tc.x[], tc.y[], tc.z[]]
    vel = [tc.v_x[], tc.v_y[], tc.v_z[]]

    pos ./= R0
    vel ./= V0
    return PhasePoint(pos, vel)
end


function to_sky(phase::PhasePoint)
    pos = phase.pos * R0
    vel = phase.vel * V0
    sc = astropy_coords.SkyCoord(x = pos[1] * u.kpc,
                                 y = pos[2] * u.kpc,
                                 z = pos[3] * u.kpc,
                                 v_x = vel[1] * u.km/u.s,
                  v_y = vel[2] * u.km/u.s,
                  v_z = vel[3] * u.km/u.s,
                  frame = galcen_frame
                 )
    tc = sc.transform_to(geocen_frame)

    return Observation(
           tc.ra[],
           tc.dec[],
           tc.distance[],
           tc.pm_ra_cosdec[],
           tc.pm_dec[],
           tc.radial_velocity[],
            )
end

function to_sky(snap::Snapshot; invert_velocity::Bool=false, verbose::Bool=false)
    observations = Observation[]

    for i in 1:length(snap)
        if verbose
            print("converting $(i)/($(length(snap))\r")
        end
        particle = snap[i]
        if invert_velocity
            vel = -1 * particle.vel
        else
            vel = particle.vel
        end
        phase = PhasePoint(particle.pos, -1 * particle.vel)
        obs = to_sky(phase)
        push!(observations, obs)
    end
    return observations
end

function rand_coord(obs::Observation, err::Observation)
    return Observation(
        ra = obs.ra,
        dec = obs.dec,
        pm_ra = obs.pm_ra + randn() * err.pm_ra,
        pm_dec = obs.pm_dec + randn() * err.pm_dec,
        radial_velocity = obs.radial_velocity + randn() * err.radial_velocity,
        distance = obs.distance + randn() * err.distance,
       )
end

function rand_coords(obs::Observation, err::Observation, N::Int)
    return [rand_coord(obs, err) for _ in 1:N]
end

function PhasePoint(x::F, y::F, z::F, vx::F, vy::F, vz::F)
    pos = [x,y,z]
    vel = [vx,vy,vz]
    return PhasePoint(pos, vel)
end

function Base.show(io::IO, pp::PhasePoint)
    @printf io "(%4.2f, %4.2f, %4.2f) kpc, " pp.pos...
    @printf io "(%4.4f, %4.4f, %4.4f) km/s" pp.vel...
    return io
end


function Base.show(io::IO, obs::Observation)
    @printf io "observation at (%4.2f, %4.2f)"  obs.ra obs.dec
    return io
end




