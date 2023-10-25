import lilguys
import numpy as np

def main():
    p0 = 1
    v0 = 0.5
    r = np.array([2, 0, 0])

    v1 = np.array([0, v0, 0])
    v2 = 0.8*v1
    v3 = 1.2*v1

    pos = np.array([r, -r, r])
    vel = np.array([v1, v2, v3])
    m = 0

    s = lilguys.Snapshot(pos, vel, m=m)
    s.save("hernquist.hdf5")


if __name__ == "__main__":
    main()
