#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 starsname"
    exit 1
fi

set -x

for f in $1/*.fits; do
    stellar_profile.jl $f -s --mass-column probability --bin-method both
done
