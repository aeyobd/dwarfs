using LilGuys
using DataFrames

snap = LilGuys.Snapshot("../orbit_smallperi/1")
N_part = length(snap)

stars_df = DataFrame(probability=1 / N_part, 
                     index = snap.index |> sort,
                    )


if !isdir("stars")
    mkdir("stars")
end


LilGuys.write_hdf5_table("stars/probabilities_stars.hdf5", stars_df, overwrite=true)
