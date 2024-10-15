#!/bin/bash

cd submodules/Move-r/
cp -rf ./patched-files/* . 

cd ../../
mkdir build
cd build
cmake ..
make
