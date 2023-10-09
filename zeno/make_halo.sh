rm -f nfw_halo nfw_1e6 nfw_1e6.txt nfw_scaled.hdf5 nfw_1e6.hdf5

$ZENOPATH/bin/halogsp nfw_halo
$ZENOPATH/bin/gspmodel nfw_halo nfw_1e6 nbody=$1
$ZENOPATH/bin/tsf nfw_1e6 maxline=100000000 > nfw_1e6.txt
python3 read_tsf.py nfw_1e6.txt nfw_1e6.hdf5
python3 process_zeno_h5.py nfw_1e6.hdf5
python rescale.py nfw_1e6.hdf5 nfw_scaled.hdf5 5 1e10
