from lilguys import Snapshot, units, profile
import argparse
import numpy as np



def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str,
            help="input hdf5 file")
    parser.add_argument('output', type=str,
            help="output hdf5 file")
    parser.add_argument('-p', '--position', type=float, nargs="+",
            help="initial position in kpc")
    parser.add_argument('-v', '--velocity', type=float, nargs="+",
            help="initial velocity in km/s")
    parser.add_argument('--max_radius', type=float, default=None,
            help="clip radius")

    args = parser.parse_args()
    return args

def set_in_orbit(snap, p, v, max_radius=None):
    centered = profile.center_snapshot(snap, verbose=True)
    if max_radius is not None:
        centered = centered.filter(centered.r < max_radius)
    dp = np.array(p) / units.R_0
    dv = np.array(v) / units.V_0
    
    centered.shift(dp, dv, inplace=True)
    pf, vf = profile.get_most_bound(centered)


    print("center at ", pf*units.R_0)
    print("moving at ", vf*units.V_0)

    return centered

def main():
    args = get_args()
        
    snap = Snapshot.file(args.input)
    centered = set_in_orbit(snap, args.position, args.velocity, args.max_radius)
    centered.save(args.output)

if __name__=="__main__":
    main()
