
rm Gadget
SCRIPT_DIR=$(pwd)
cd ~/gadget4 || exit

make CONFIG=$SCRIPT_DIR/Config.sh BUILD_DIR=$SCRIPT_DIR/build EXEC=$SCRIPT_DIR/Gadget
