setenv MC_WD $PWD
cd src
make
mv membrane.exe $MC_WD
cd $MC_WD
if ( -f gen_input.py ) then
  echo "You have  a copy of gen_input.py."
else
  echo " Copying python file to generate input files in this directory."
  cp $MC_HOME/bin/gen_input.py .
  echo "---- DONE ----"
endif
echo "Create input files using this python script."
echo "--------------------------------------------"
