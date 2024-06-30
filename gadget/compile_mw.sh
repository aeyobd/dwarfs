
rm GadgetMW
SCRIPT_DIR=$(pwd)
cd $GADGET_SOURCE || exit

make CONFIG=$SCRIPT_DIR/ConfigMW.sh BUILD_DIR=$SCRIPT_DIR/build_mw EXEC=$SCRIPT_DIR/GadgetMW
