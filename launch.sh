#!/bin/bash

alpha=8
fileLength=200000

./build/build_MLF example $alpha

for index in 8 16 32 64; do
    ./build/MOVE_EDSBWTSearch example ../datasets/randomPatterns/$fileLength"_"$index.txt
done