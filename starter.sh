#!/bin/bash

set -e   # exit on any command failure

INCLUDE_PATH="-I/usr/bin/eigen"
LIBS="-lfmt"

sources=("src/lodepng.cpp" "src/utils.cpp" "src/powerIteration.cpp" "src/main.cpp")
executables=()

for src in "${sources[@]}"; do
    exe="${src%.cpp}.o"
    g++ $INCLUDE_PATH -c "$src" $LIBS -o "$exe"
    executables+=("$exe")
done

g++ "${executables[@]}" $LIBS -o main.exe
echo "main.cpp linked and compiled successfully"
./main.exe