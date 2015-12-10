#!/bin/bash

SOURCE=../trunk
CURRENT=$(pwd)

rm *windows*.zip

cd ${SOURCE}
make clean
#make distclean
make compiler=mingw64
cd ${CURRENT}

. build-mingw64.sh
