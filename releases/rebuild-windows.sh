#!/bin/bash

SOURCE=../trunk
CURRENT=$(pwd)

rm *windows*.zip

cd ${SOURCE}
make clean
#make distclean
make build=static
cd ${CURRENT}

. build-windows.sh
