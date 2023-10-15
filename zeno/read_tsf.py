import h5py
import sys
import re
from os import path
import numpy as np

from lilguys import Snapshot


def read_file(filename):
    header = True
    table = {}
    sublist = []

    variable = None
    ndims = None

    with open(filename) as file:
        for line in file:
            lintems = line.strip().split(" ")
            if header:
                if len(lintems) < 2:
                    continue

                if "Particles" in lintems[1]:
                    header = False
                    print(lintems)
                    continue
                continue


            match = re.search(r"[a-zA-Z]+", lintems[0])
            if match:
                kind = match.group()
                if kind == "float":
                    variable = re.search(r"[a-zA-Z]+", lintems[1]).group()
                    dims = re.findall(r"\d+", lintems[1])
                    dims = [float(d) for d in dims]
                    if len(dims) == 1:
                        ndims = 1
                    elif len(dims) == 2:
                        ndims = dims[1]
                    else:
                        print(dims)
                        raise Exception("only works for 2d data")

                    N = dims[0]
                    print(variable)
                    print(N, ndims)

                    table[variable] = [float(item) for item in lintems[2:]]

                elif kind == "tes":
                    break
                elif kind == "e": #scientific exponent is numeric
                    for item in lintems:
                        table[variable].append(float(item))
                else:
                    raise Exception("unknown type: ", kind)
            else:
                for item in lintems:
                    table[variable].append(float(item))


    return table


def inflate_list(l, ndim):
    if len(l) % ndim != 0:
        raise Exception("dimension mismatch: ", len(l), ndim)
    N = len(l)//ndim
    out = []
    for i in range(N):
        lix = i*ndim
        uix = (i+1)*ndim 
        out.append(l[lix:uix])
    return np.array(out)


def make_header(N, mass):
    Npart = 2

    header = {}
    header["NumPart_ThisFile"] = (0, N)
    header["NumPart_Total"] = (0, N)
    header["MassTable"] = (0., mass)
    header["Time"] = 0.
    header["Redshift"] = 0.
    header["BoxSize"] = 350.
    header["NumFilesPerSnapshot"] = 1
    return header

def output_to_hdf5(table, filename):
    N = len(table["Mass"])
    mass = table["Mass"][0]
    print(mass)
    header = make_header(N, mass)
    snap = Snapshot(
            table["Position"],
            table["Velocity"],
            IDs = np.arange(N),
            header=header
            )
    # snap.gadget4 = True
    snap.save(filename)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("requires filename, outname")
        sys.exit(1)
    filename = sys.argv[1]
    outname = sys.argv[2]

    particles = read_file(filename)

    particles["Position"] = inflate_list(particles["Position"], 3)
    particles["Velocity"] = inflate_list(particles["Velocity"], 3)
    for var, l in particles.items():
        print(var, len(l),  l[:5])
    output_to_hdf5(particles, outname)
