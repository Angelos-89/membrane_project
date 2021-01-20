#!/bin/bash
export MC_WD=PWD        
mkdir -p src
cd src
ln -s /home/angelos/GIT/membrane_project/src/*.cpp .
ln -s /home/angelos/GIT/membrane_project/src/*.hpp .
ln -s /home/angelos/GIT/membrane_project/src/makefile .
ln -s /home/angelos/GIT/membrane_project/config/config.angelos .
