import h5py
import sys
import re
from os import path
import numpy as np


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
    return out

def output_to_hdf5(table, filename):
    N = len(table["Mass"])
    with h5py.File(filename, "w") as f:
        grp = f.create_group("PartType1")
        header = f.create_group("Header")
        params = f.create_group("Parameters")
        conf = f.create_group("Config")
        # m = f.create_dataset("Mass", (N,), dtype="f")
        # m[:] = table["Mass"]
        x = grp.create_dataset("Position", (N,3), dtype="f")
        x[:,:] = table["Position"]
        v = grp.create_dataset("Velocity", (N,3), dtype="f")
        v[:,:] = table["Velocity"]
        idx = grp.create_dataset("ParticleIDs", (N), dtype="i")
        idx[:] = np.arange(N)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("requires filename")
        sys.exit(1)
    filename = sys.argv[1]

    particles = read_file(filename)

    particles["Position"] = inflate_list(particles["Position"], 3)
    particles["Velocity"] = inflate_list(particles["Velocity"], 3)
    for var, l in particles.items():
        print(var, len(l),  l[:5])
    hdfname = path.splitext(filename)[0]
    hdfname += ".hdf5"
    output_to_hdf5(particles, hdfname)
