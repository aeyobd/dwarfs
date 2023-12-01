from lilguys import Snapshot, units
import argparse
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('output')
    parser.add_argument('-M', '--mass', type=float)
    parser.add_argument('-R', '--radius', type=float)
    parser.add_argument('--max_radius', type=float, default=None)

    args = parser.parse_args()
    return args


def main():
    args = get_args()
        
    snap = Snapshot.file(args.input)
    r_scale = args.radius / units.R_0
    m_scale = args.mass / units.M_0
    v_scale = np.sqrt(units.G * m_scale / r_scale)
    scaled = snap.scale(r_scale, v_scale, m_scale)
    if args.max_radius:
        scaled = scaled.filter(scaled.r < args.max_radius)
    # scaled.gadget4 = True
    scaled.save(args.output)

if __name__=="__main__":
    main()
