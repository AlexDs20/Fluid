#! /bin/bash

CONFIG="${1:-release}"

if [ -d "build/bin/" ]; then
    rm -r build/bin/{Fluid,Tests}
fi

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

FAILED=0

if [ $SUCCESS -eq 0 ]; then
    for i in {1..1}
    do
        echo "Run: $i"
        timeout 2s ./build/bin/Fluid/$CONFIG/Fluid
        if [ $? -ne 124 ]; then
            ((FAILED++))
        fi
        echo "Failed runs: $FAILED / $i"
    done
else
    echo "Build failed"
fi
