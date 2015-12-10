#!/bin/bash

SOURCE=../trunk
SPS=./sps/3.0/linux_x64
SPECPLOT=./specplot/2.0/linux_x64
CURRENT=$(pwd)

echo ${CURRENT}

# sps
cd ${SPS}
rm -rf *.zip
. build-static.sh
. build-dynamic.sh
cp *.zip ${CURRENT}
cd ${CURRENT}


# specplot
cd ${SPECPLOT}
rm -rf *.zip
. build-static.sh
. build-dynamic.sh
cp *.zip ${CURRENT}
cd ${CURRENT}

