setenv MC_WD $PWD
cd src
make r_spec.exe 
mv r_spec.exe $MC_WD
cd $MC_WD
echo "------DONE now test fft-----------"
