#!/bin/sh

g++ -I/$HOME/.local/include/ -std=c++11 -O3 -DNDEBUG "test.cpp" -o "test.exe"
./"test.exe"
