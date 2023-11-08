import lilguys
import numpy as np

def main():
    r0 = 20
    v0 = 0.15
    r = np.array([r0, 0, 0])

    v1 = np.array([0, v0, 0])
    v2 = 0.5*v1
    v3 = 1.2*v1

    pos = np.array([r, -r, r])
    vel = np.array([v1, v2, v3])
    m = 0

    s = lilguys.Snapshot(pos, vel, m=m)
    s.save("nfw_test.hdf5")


if __name__ == "__main__":
    main()
