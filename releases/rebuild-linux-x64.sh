#!/bin/bash

SOURCE=../trunk
CURRENT=$(pwd)

rm *-x64.zip

cd ${SOURCE}
make clean
#make distclean
make build=static
cd ${CURRENT}

. build-linux-x64.sh
