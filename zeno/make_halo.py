import os
import argparse
from read_tsf import tsf_to_hdf5


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('N', type=str)

    args = parser.parse_args()
    return args

def main():
    args = get_args()
    N = int(float(args.N))
    N_name = args.N

    print(f"using {N_name} particles")

    print("cleaning up")
    halo_path = 'nfw_halo.gsp'
    txt_path = f"nfw_{N_name}.txt"
    model_path = f'nfw_{N_name}.gsp'
    hdf5_path = f'nfw_{N_name}.hdf5'
    os.system(f"rm {halo_path} {model_path} {txt_path} {hdf5_path}")

    print('generating profile ')
    os.system(f'$ZENOPATH/bin/halogsp {halo_path}')

    print('generating nbody model')
    os.system(f'$ZENOPATH/bin/gspmodel {halo_path} {model_path} nbody={N}')


    print("writing to text file")
    os.system(f"$ZENOPATH/bin/tsf {model_path} maxline={N} > {txt_path}")


    print("converting to hdf5")
    tsf_to_hdf5(txt_path, hdf5_path)


if __name__ == "__main__":
    main()

