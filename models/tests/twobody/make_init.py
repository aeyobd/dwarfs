import lilguys
import numpy as np

def main():
    p0 = 1
    v0 = 0.3
    x0 = np.array([[p0,0,0], 
        [-p0,0,0]])
    v0 = np.array([[0,v0,0], 
        [0,-v0,0]])

    m = 1

    s = lilguys.Snapshot(x0, v0, m=m)
    s.save("twobody.hdf5")


if __name__ == "__main__":
    main()
