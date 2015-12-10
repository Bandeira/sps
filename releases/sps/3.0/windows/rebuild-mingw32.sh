rm -rf sps
rm *.zip

cd ../../../../trunk
make clean
make compiler=mingw32


cd ../releases/sps/3.0/windows

. build-mingw32.sh

