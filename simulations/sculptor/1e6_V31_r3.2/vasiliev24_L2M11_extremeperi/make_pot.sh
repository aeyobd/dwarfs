
pot_dir=$DWARFS_ROOT/agama/potentials/vasiliev24/L2M11/

rm agama_potential.ini
rm trajlmc.txt boundmass.txt accel.txt
rm -rf pot

cp $pot_dir/potential.ini agama_potential.ini
cp $pot_dir/trajlmc.txt .
cp $pot_dir/boundmass.txt .
cp $pot_dir/accel.txt .
cp -r $pot_dir/pot .

