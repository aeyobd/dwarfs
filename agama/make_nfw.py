import agama
from astropy.io import fits
import argparse
from math import pi


def main():
    args = parse_args()
    pos, vel, mass = sample_nfw(args.N)
    write_to_fits(pos, vel, mass, args.output)

def parse_args():
    parser = argparse.ArgumentParser(description='Create NFW potential')
    parser.add_argument('-n', '--N', type=float, default=1e4, help='Number of particles')
    parser.add_argument('-o', '--output', type=str, default='nfw.fits', help='Output file name')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('-t', '--truncation_radius', type=float, default=16, help='Truncation radius over scale radius')

    return parser.parse_args()



def sample_nfw(N=1e4):
    # Create a NFW potential
    if not isinstance(N, (int, float)):
        raise ValueError('N must be a number')


    if N % 1 != 0:
        raise ValueError('N must be an integer')

    if N < 1:
        raise ValueError('N must be greater than 1')

    N = int(N)
    scaleRadius = 1
    mass = 1
    cutoff = 16

    rho0 = mass / (4*pi * scaleRadius**3 ) / 3

    pot = agama.Potential(type='Spheroid', densityNorm=rho0,
        gamma=1, beta=3, alpha=1, scaleRadius=scaleRadius,
        outerCutoffRadius = cutoff*scaleRadius,
        cutoffStrength=1)

    df = agama.DistributionFunction(type="QuasiSpherical", potential=pot)
    gm = agama.GalaxyModel(pot, df)

    posvel, mass = gm.sample(N)
    pos = posvel[:,:3]
    vel = posvel[:,3:]

    return pos, vel, mass



def write_to_fits(pos, vel, mass, output='nfw.fits'):
    hdu = fits.HDUList()
    hdu.append(fits.PrimaryHDU())
    hdu[0].header['N'] = len(pos)
    hdu.append(fits.BinTableHDU.from_columns([
        fits.Column(name='x', format='D', array=pos[:,0]),
        fits.Column(name='y', format='D', array=pos[:,1]),
        fits.Column(name='z', format='D', array=pos[:,2]),
        fits.Column(name='vx', format='D', array=vel[:,0]),
        fits.Column(name='vy', format='D', array=vel[:,1]),
        fits.Column(name='vz', format='D', array=vel[:,2]),
        fits.Column(name='mass', format='D', array=mass),
    ]))
    hdu.writeto(output, overwrite=True)
    print(f'Wrote {len(pos)} particles to {output}')



if __name__ == '__main__':
    main()

