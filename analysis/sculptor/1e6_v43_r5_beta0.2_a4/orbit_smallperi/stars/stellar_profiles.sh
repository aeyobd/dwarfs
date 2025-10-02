#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 starsname"
    exit 1
fi

for f in $1/*.fits; do
    stellar_profile.jl $f --mass-column probability --bin-method both --centre first
done
