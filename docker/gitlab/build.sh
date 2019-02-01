# TODO remove net line
env
module restore dtk
# FIXME CI_PROJECT_DIR env variable was not set
DTK_DIR=$PWD
ln -sf $DTK_DIR $TRILINOS_DIR
mkdir build && cd build
$DTK_DIR/scripts/olcf_ascent_cmake
source $DTK_DIR/scripts/set_kokkos_env.sh
make -j36
