#! /bin/bash

CONFIG=$1

if [ ! -d "generated" ]; then
    premake5 gmake2
fi

pushd "generated"
    #make clean
    make -j$(nproc) config=$CONFIG
    SUCCESS=$?
popd

if [ $SUCCESS -eq 0 ]; then
    ./build/bin/Fluid/$CONFIG/Fluid
fi
