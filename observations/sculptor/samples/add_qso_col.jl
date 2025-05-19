using PyFITS
using DataFrames


df = read_fits("../data/jensen+24_wide_2c.fits")
df_gaia = read_fits("../data/gaia_4deg_cen.fits")[:, [:source_id, :in_qso_candidates, :in_galaxy_candidates]]


leftjoin!(df, df_gaia, on=:source_id)


@info names(df)[[1, end-2, end-1, end]]
write_fits("jensen+24_wide_2c_extracols.fits", df)
