import lilguys
import numpy as np

def main():
    p0 = 1
    v0 = 0.3
    r1 = np.array([0.97000436, -0.24308753, 0])
    r2 = -r1
    r3 = np.zeros(3)

    v3 = np.array([-0.93240737, -0.86473146, 0])
    v1 = -v3/2
    v2 = -v3/2

    pos = np.array([r1, r2, r3])
    vel = np.array([v1, v2, v3])
    m = 1

    s = lilguys.Snapshot(pos, vel, m=m)
    s.save("threebody.hdf5")


if __name__ == "__main__":
    main()
