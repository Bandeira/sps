rm -rf sps
rm *.zip

cd ../../../../trunk
make clean
make build=static


cd ../releases/sps/3.0/linux_x64

. build-static.sh
. build-dynamic.sh


