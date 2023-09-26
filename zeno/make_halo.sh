rm -f nfw_halo nfw_1e6 nfw_1e6.txt

halogsp nfw_halo
gspmodel nfw_halo nfw_1e6 nbody=1000000
tsf nfw_1e6 maxline=100000000 > nfw_1e6.txt
python read_tsf.py nfw_1e6.txt
