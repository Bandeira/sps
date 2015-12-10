rm -rf sps
rm *.zip

cd ../../../../trunk
make clean
make build=static


cd ../releases/sps/3.0/windows

. build.sh

