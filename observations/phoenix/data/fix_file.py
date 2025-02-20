from astropy.table import Table
import numpy as np

t = Table.read("j24_1c_err.fits")

cols = t.colnames
cols = [col for col in cols if col != "IDVEL"]


t[cols].write("j24_1c.fits", overwrite=True)
