#! /bin/bash

CONFIG="${1:-release}"

rm -r build/bin/{Fluid,Tests}

if [ "premake5.lua" -nt "generated/Makefile" ]; then
    echo "Newer premake, deleting generated"
    rm -r "generated"
fi

if [ ! -d "generated" ]; then
    echo "premake5 gmake2"
    premake5 gmake2
fi

pushd "generated"
    make -j$(nproc) config=$CONFIG
    SUCCESS=$?
popd

if [ $SUCCESS -eq 0 ]; then
    ./build/bin/Fluid/$CONFIG/Fluid
else
    echo "Build failed"
fi
