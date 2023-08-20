#! /bin/bash

CONFIG=$1

rm -rf build generated

premake5 gmake2

pushd generated
make -j$(nproc) config=$CONFIG
popd

./build/bin/Fluid/$CONFIG/Fluid
