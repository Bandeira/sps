#!/bin/bash

SOURCE=../trunk
CURRENT=$(pwd)

rm *-x86.zip

cd ${SOURCE}
make clean
#make distclean
make build=static type=32
cd ${CURRENT}

. build-linux-x86.sh
