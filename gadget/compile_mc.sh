rm GadgetMC

SCRIPT_DIR=$(pwd)
cd $GADGET_SOURCE || exit

make CONFIG=$SCRIPT_DIR/ConfigMC.sh BUILD_DIR=$SCRIPT_DIR/build_mc EXEC=$SCRIPT_DIR/GadgetMC
