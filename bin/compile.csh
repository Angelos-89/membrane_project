setenv MC_WD $PWD
cd src
make
mv membrane.exe $MC_WD
cd $MC_WD
if ( -f gen_input.py ) then
  echo "you have  a copy of gen_input"
else
  echo " copying python file to generate input file in theis directory"
  cp $MC_HOME/bin/gen_input.py .
  echo "------done ----"
endif
echo "Create input files using this python script"
echo "------------------------------------"
