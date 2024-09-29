#!/bin/bash
#Updates file access times for a particular file in /scratch/t/todelete/current/
#Usage: sh touchfiles.sh <file location>
#Example: sh touchfiles.sh /scratch/t/todelete/current/3002243__fherwig_____fherwig________6.01T___626718files

counter=0
while read -r line; do
    linearr=($line)
    echo -en "\r File number ${counter}; ${linearr[12]}"
    touch "${linearr[12]}"
    counter=$(( $counter+1 ))
done < "$1"
