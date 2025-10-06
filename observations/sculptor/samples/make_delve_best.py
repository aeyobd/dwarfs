from astropy.table import Table


df = Table.read("../data/delve_best.fits")

filt = df["flags_g"] <= 4
filt &= df["flags_r"] <= 4
filt &= df["extended_class_g"] <= 1
filt &= df["extended_class_g"] >= 0

df_filtered = df[filt]


df_filtered.write("delve_best_sample.fits")
