#!/bin/bash

SOURCE=../trunk
SPS=./sps/3.0/windows
SPECPLOT=./specplot/2.0/windows
CURRENT=$(pwd)

echo ${CURRENT}

# sps
cd ${SPS}
rm -rf *.zip
. build.sh
cp *.zip ${CURRENT}
cd ${CURRENT}


# specplot
cd ${SPECPLOT}
rm -rf *.zip
. build.sh
cp *.zip ${CURRENT}
cd ${CURRENT}

