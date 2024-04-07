#!/bin/bash

# read in argument
if [ $# -ne 1 ]; then
  echo "Usage: $0 <N>"
  exit 1
fi

set -e

printf -v N '%0.0f' $1
filename="nfw_$1"

echo using $N particles

halo_path=nfw_halo.gsp
model_path=$filename.gsp
txt_path=$filename.txt
hdf5_path=${filename}_uncentred.hdf5
out_path=$filename.hdf5


echo cleaning up

rm -f $halo_path $model_path $txt_path $hdf5_path $out_path

echo generating halo
$ZENOPATH/bin/halogsp $halo_path

echo generating nbody model
$ZENOPATH/bin/gspmodel $halo_path $model_path nbody=$N

echo writing to text file
$ZENOPATH/bin/tsf $model_path maxline=$N > $txt_path

echo converting to hdf5

julia parse_zeno.jl $txt_path $hdf5_path

echo centreing

../scripts/centre_snapshot.jl -c $hdf5_path $out_path
echo completed
