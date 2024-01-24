# rescales input relative to fornax
python3 ~/dwarfs/set_in_orbit.py isolation.hdf5 orbit_1.hdf5 -p 31.57 102.59 122.76 -v -37.8505 105.2257 -78.1885 --max_radius 50
python3 ~/dwarfs/set_in_orbit.py isolation.hdf5 orbit_2.hdf5 -p -44.02 -7.84 -135.08 -v 7.6355 -102.3527 3.4242 --max_radius 50
python3 ~/dwarfs/set_in_orbit.py isolation.hdf5 orbit_3.hdf5 -p 42.69 -167.29 38.24 -v 41.7971 40.7130 131.5140 --max_radius 50
