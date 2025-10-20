#!/bin/bash

filename=$1
outname="${filename%.*}.table"

rm $outname
julia $filename > $outname

