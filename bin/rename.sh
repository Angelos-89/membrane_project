#---------------------------------------#
#This script is used to rename the files hfield_0.h5 to hfield_eq_0.h5
#so that we can start a simulation from an earlier one.
#Uncomment the following code and change the directory according
#to the .h5 files location.
#---------------------------------------#

#FILEPATH="/home/angelos-89/bash_scripting"
#SUBSTR2="eq_"
#SUBSTR4=".h5"
#for FILE in $FILEPATH/*.h5
#do
#    FILENAME=${FILE##*/}
#    VAR1=${FILENAME#*_}
#    SUBSTR3=${VAR1%???}
#    SUBSTR1="${FILENAME:0:7}"
#    mv $FILENAME $SUBSTR1$SUBSTR2$SUBSTR3$SUBSTR4
#done
