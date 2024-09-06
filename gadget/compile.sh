#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 <config_name>"
    exit 1
fi

paramname=$1
prefix="Config"
suffix=".sh"
name=${paramname#"$prefix"}
name=${name%"$suffix"}

rm "Gadget$name"
SCRIPT_DIR=$(pwd)
cd $GADGET_SOURCE || exit

make CONFIG=$SCRIPT_DIR/$paramname BUILD_DIR="$SCRIPT_DIR/build_$name" EXEC="$SCRIPT_DIR/Gadget$name"
