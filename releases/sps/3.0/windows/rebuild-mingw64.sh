rm -rf sps
rm *.zip

cd ../../../../trunk
make clean
make compiler=mingw64


cd ../releases/sps/3.0/windows

. build-mingw64.sh

