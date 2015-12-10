rm -rf specplot
rm *.zip

cd ../../../../trunk
make clean
make compiler=mingw32


cd ../releases/specplot/2.0/windows

. build-mingw32.sh

