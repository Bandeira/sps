rm -rf sps
rm *.zip

cd ../../../../trunk
make clean
make type=32 build=static

cd ../releases/sps/3.0/linux_x86

. build-static.sh
. build-dynamic.sh

