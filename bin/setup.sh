#!/bin/bash
export MC_WD=PWD        
mkdir -p src
cd src
ln -s /home/angelos/GIT/MC_membrane/src/*.cpp .
ln -s /home/angelos/GIT/MC_membrane/src/*.h .
ln -s /home/angelos/GIT/MC_membrane/src/makefile.der_test .
ln -s /home/angelos/GIT/MC_membrane/config/config.angelos .
