#! /bin/bash

CONFIG=$1

rm -r build/bin/Fluid/

if [ ! -d "generated" ]; then
    premake5 gmake2
fi

pushd "generated"
    make -j$(nproc) config=$CONFIG
    SUCCESS=$?
popd

if [ $SUCCESS -eq 0 ]; then
    ./build/bin/Fluid/$CONFIG/Fluid
fi
