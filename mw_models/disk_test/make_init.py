import lilguys
import numpy as np

def main():
    p0 = 1
    v0 = 0.8
    r = np.array([8, 0, 0])

    v1 = np.array([0, v0, 0.05*v0])
    v2 = np.array([0, 1.2*v0, 0.05*v0])
    v3 = np.array([0, 0.8*v0, 0.05*v0])

    pos = np.array([r, -r, r])
    vel = np.array([v1, v2, v3])
    m = 0

    s = lilguys.Snapshot(pos, vel, m=m)
    s.save("disk_test.hdf5")


if __name__ == "__main__":
    main()
