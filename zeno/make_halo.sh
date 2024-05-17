#!/bin/bash

# read in argument
if [ $# -ne 3 ]; then
  echo "Usage: $0 <N> <halo.gsp> <out>"
  exit 1
fi

set -e

printf -v N '%0.0f' $1
halo_path=$2
filename=$3

echo using $N particles

model_path=$filename.gsp
txt_path=$filename.txt
hdf5_path=${filename}_uncentred.hdf5
out_path=$filename.hdf5


echo cleaning up

rm -f  $model_path $txt_path $hdf5_path $out_path

echo generating nbody model
echo reading $halo_path as profile table
echo writing to $model_path
$ZENOPATH/bin/gspmodel $halo_path $model_path nbody=$N

echo writing $model_path to text file $txt_path
$ZENOPATH/bin/tsf $model_path maxline=$N > $txt_path

echo converting $txt_path to hdf5 file $hdf5_path

julia parse_zeno.jl $txt_path $hdf5_path

echo centreing $hdf5_path and writing to $out_path
../scripts/centre_snapshot.jl $hdf5_path $out_path -m com 
echo completed
