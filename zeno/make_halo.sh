rm -f nfw_halo.gsp nfw_1e4.gsp nfw_1e4.dat nfw_scaled.hdf5 nfw_1e4.hdf5

$ZENOPATH/bin/halogsp nfw_halo.gsp
$ZENOPATH/bin/gspmodel nfw_halo.gsp nfw_1e4.gsp nbody=10000
$ZENOPATH/bin/tsf nfw_1e4.gsp maxline=100000000 > nfw_1e4.dat
python3 read_tsf.py nfw_1e4.dat nfw_1e4.hdf5
# python3 process_zeno_h5.py nfw_1e4.hdf5
# python3 rescale.py nfw_1e4.hdf5 nfw_scaled.hdf5 2 2e8
# python3 cent.py nfw_1e4.hdf5 nfw_scaled.hdf5 2 2e8
