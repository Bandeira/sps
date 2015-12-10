rm -rf specplot
rm *.zip

cd ../../../../trunk
make clean
make build=static


cd ../releases/specplot/2.0/windows

. build.sh

