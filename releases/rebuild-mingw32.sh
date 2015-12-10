#!/bin/bash

SOURCE=../trunk
CURRENT=$(pwd)

rm *windows*.zip

cd ${SOURCE}
make clean
#make distclean
make compiler=mingw32
cd ${CURRENT}

. build-mingw32.sh
