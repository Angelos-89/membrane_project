#!/bin/bash
export MC_WD=PWD        
mkdir -p src
cd src
ln -s $MC_HOME/src/*.cpp .
ln -s $MC_HOME/src/*.hpp .
ln -s $MC_HOME/src/makefile .
