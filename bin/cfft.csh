setenv MC_WD $PWD
cd src
make ftest
mv ftest.exe $MC_WD
cd $MC_WD
echo "------DONE now run-----------"
