#!/bin/bash

cd Move-r/
cp -rf ./patched-files/* . 

cd ../
mkdir build
cd build
cmake ..
make
