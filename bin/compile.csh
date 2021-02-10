setenv MC_WD $PWD
cd src
make
mv membrane.exe $MC_WD
cd $MC_WD
echo "------------------------------------"
echo "Make sure you have an input file in this directory"
echo "------------------------------------"
