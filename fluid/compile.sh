#! /bin/bash

CONFIG=$1

premake5 gmake2

pushd generated
make clean
make -j$(nproc) config=$CONFIG
popd

./build/bin/Fluid/$CONFIG/Fluid
