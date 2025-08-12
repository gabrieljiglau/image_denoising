#!/bin/bash

set -e   # exit on any command failure

lode_png=src/lodepng.cpp
lode_png_out=src/lodepng.o
if [ ! -f "$lode_png_out" ]; then
   g++ -c "$lode_png" -o "$lode_png_out"
   if [ $? -ne 0 ]; then
      echo "Compiling "$lode_png" failed"
      exit 1
   fi
fi

utils_path=src/utils.cpp
utils_out=src/utils.o
g++ -c "$utils_path" -lfmt -o "$utils_out"
if [ $? -ne 0 ]; then
   echo "Compiling $utils_path failed"
   exit 2
fi


power_iteration_path=src/powerIteration.cpp
power_iteration_out=src/powerIteration.o
g++ "-I/usr/bin/eigen" -c "$power_iteration_path" -lfmt -o "$power_iteration_out"
if [ $? -ne 0 ]; then
   echo "Compiling $power_iteration_path failed"
   exit 3
fi


main=src/main.cpp
main_out=src/main.o
g++ "-I/usr/bin/eigen" -c $main -o "$main_out" 
if [ $? -ne 0 ]; then
   echo "Compiling $main failed"  
   exit 4
fi

g++ "$main_out" "$lode_png_out" "$utils_out" "$power_iteration_out" -lfmt -o main.exe
echo "main.cpp linked and compiled successfully"
./main.exe