rm -rf specplot
rm *.zip

cd ../../../../trunk
make clean
make compiler=mingw64


cd ../releases/specplot/2.0/windows

. build-mingw64.sh

