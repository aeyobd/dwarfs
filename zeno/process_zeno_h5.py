import h5py
import sys



def process(filename):

    with h5py.File(filename, "a") as f:
        header = f["Header"]
        attr = header.attrs

        l = len(f["PartType1"]["ParticleIDs"])
        M = 0.3
        Npart = 2

        attr.create("NumPart_ThisFile", [0,l], (Npart), dtype="uint")
        attr.create("NumPart_Total", [0,l], (Npart), dtype="uint")

        attr.create("MassTable", [0, M], (Npart), dtype="double")
        attr.create("Time", 0, (1), dtype="double")
        attr.create("Redshift", 0, (1), dtype="double")
        attr.create("BoxSize", 0, (1), dtype="double")
        attr.create("NumFilesPerSnapshot", (1), dtype="int")


if __name__ == "__main__":
    process(sys.argv[1])
