
rm GadgetMW
SCRIPT_DIR=$(pwd)
cd ~/gadget4 || exit

make CONFIG=$SCRIPT_DIR/ConfigMW.sh BUILD_DIR=$SCRIPT_DIR/build_mw EXEC=$SCRIPT_DIR/GadgetMW
