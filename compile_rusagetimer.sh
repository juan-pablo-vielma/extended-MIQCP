#!/bin/sh
gcc -c -fPIC rusagetimer.c -o rusagetimer.o
gcc -shared -W1,-soname,rusagetimer.so.1 -o rusagetimer.so.1.0.1 rusagetimer.o
g++ -dynamiclib -undefined suppress -flat_namespace rusagetimer.o -o rusagetimer.dylib
