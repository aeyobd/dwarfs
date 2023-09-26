import h5py
import sys





def process(filename):
    with h5py.File(filename, "a") as f:
        header = f["Header"]


if __name__ == "__main__":
    process(sys.argv[1])
