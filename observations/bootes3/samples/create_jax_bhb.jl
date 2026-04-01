using PyFITS

stars = read_fits("../data/j24_1c.fits")
stars[:, "L_CMD_SAT"] = stars.L_CMD_SAT_HB
filt = stars.L_CMD_SAT ./ stars.L_CMD_BKD .> 0
filt .&= stars.F_BEST .== 1

stars_out = stars[filt, :]

stars_out[:, :L_sat] = @. stars_out.L_CMD_SAT * stars_out.L_PM_SAT
stars_out[:, :L_bg] = @. stars_out.L_CMD_BKD * stars_out.L_PM_BKD
write_fits("j24_bhb.fits", stars_out, overwrite=true)
