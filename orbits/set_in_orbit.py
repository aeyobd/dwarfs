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

    args = parser.parse_args()
    return args


def main():
    args = get_args()
        
    snap = Snapshot.file(args.input)

    centered = profile.center_snapshot(snap, verbose=True)
    dp = np.array(args.position) / units.R_0
    dv = np.array(args.velocity) / units.V_0
    
    centered.shift(dp, dv, inplace=True)
    pf, vf = profile.get_most_bound(centered)
    print("center at ", pf*units.R_0)
    print("moving at ", vf*units.V_0)

    centered.save(args.output)

if __name__=="__main__":
    main()
