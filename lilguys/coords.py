import numpy as np
from astropy.coordinates import Galactocentric, ICRS, SkyCoord
from astropy import units as u
from . import units


galcen_frame = Galactocentric(
    galcen_distance = 8.29*u.kpc,
    z_sun=0*u.pc,
    galcen_v_sun = [11.1, 240.3+12.24, 7.25] * (u.km/u.s)
)

geocen_frame = ICRS()


class observation:
    def __init__(self, ra, dec, 
            distance=None, pm_ra=None, pm_dec=None, 
            radial_velocity=None,
            distance_err=None, pm_ra_err=None, pm_dec_err=None,
            radial_velocity_err=None):
        self.ra = ra
        self.dec = dec
        self.distance = distance
        self.pm_ra = pm_ra
        self.pm_dec = pm_dec
        self.radial_velocity = radial_velocity
        self.distance_err = distance_err
        self.pm_ra_err = pm_ra_err
        self.pm_dec_err = pm_dec_err
        self.radial_velocity_err = radial_velocity_err

    def __str__(self):
        return f"observation at ({self.ra}, {self.dec})"

    def __repr__(self):
        return str(self)


class phase_point:
    def __init__(self, x, y, z, v_x, v_y, v_z):
        self.x = x
        self.y = y 
        self.z = z
        self.v_x = v_x
        self.v_y = v_y
        self.v_z = v_z
    
    @property
    def v(self):
        return [self.v_x, self.v_y, self.v_z]
    
    @property 
    def pos(self):
        return [self.x, self.y, self.z]
    
    def __str__(self):
        s = "["
        s += f"{self.x:0.2f}, "
        s += f"{self.y:0.2f}, "
        s += f"{self.z:0.2f}"
        s += "]kpc, \t["
        s += f"{self.v_x:0.4f}, "
        s += f"{self.v_y:0.4f}, "
        s += f"{self.v_z:0.4f}"
        s += "] km/s"
        return s
        
    def __repr__(self):
        return str(self)


def to_galcen(obs):
    sc = SkyCoord(ra = obs.ra * u.deg,
                  dec=obs.dec * u.deg,
                  distance = obs.distance * u.kpc,
                  radial_velocity = obs.radial_velocity * u.km/u.s,
                  pm_ra_cosdec = obs.pm_ra * u.mas/u.yr,
                  pm_dec = obs.pm_dec * u.mas/u.yr
                 )
    tc = sc.transform_to(galcen_frame)
    x, y, z = [tc.x, tc.y, tc.z]
    vx, vy, vz = [tc.v_x, tc.v_y, tc.v_z]
    x = x.to("kpc").value
    y = y.to("kpc").value
    z = z.to("kpc").value
 
    vx = vx.to("km/s").value
    vy = vy.to("km/s").value
    vz = vz.to("km/s").value
    return phase_point(x, y, z, vx, vy, vz)



def to_sky(phase):
    sc = SkyCoord(
            x = phase.x * u.kpc, 
            y = phase.y * u.kpc,
            z = phase.z * u.kpc,
            v_x = phase.v_x * u.km/u.s,
            v_y = phase.v_y * u.km/u.s,
            v_z = phase.v_z * u.km/u.s,
            frame = galcen_frame
            )
    tc = sc.transform_to(geocen_frame)

    return observation(
            tc.ra.to("degree").value,
            tc.dec.to("degree").value,
            tc.distance.to("kpc").value,
            tc.pm_ra_cosdec.to("mas/yr").value,
            tc.pm_dec.to("mas/yr").value,
            tc.radial_velocity.to("km/s").value,
            )


def rand_coords(obs, N):
    """Draws N random coordinates from a given observation obs"""
    ra = np.full(N, obs.ra)
    dec = np.full(N, obs.dec)
    distance = np.random.normal(obs.distance,obs.distance_err, N)
    pm_ra = np.random.normal(obs.pm_ra,obs.pm_ra_err, N)
    pm_dec = np.random.normal(obs.pm_dec,obs.pm_dec_err, N)
    radial_velocity = np.random.normal(obs.radial_velocity,obs.radial_velocity_err, N)
    X = np.array([ra, dec, distance, pm_ra, pm_dec, radial_velocity])
    rand_obs = [observation(*x) for x in X.transpose()]
    return rand_obs
