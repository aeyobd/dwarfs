using PyFITS

stars = read_fits("../data/delve_6deg_good.fits")
filt = stars.gmag .- stars.rmag .< 0.2
filt .&= 18 .< stars.gmag .< 19.5

write_fits("../data/delve_bhb.fits", stars_out, overwrite=true)
