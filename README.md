# L-SRTDE algorithm for CEC 2024 benchmark

C++ implementation of the L-SRTDE algorithm (Linear population size reduction Success RaTe-based adaptive Differential Evolution) for the Congress on Evolutionary Computation competition on single-objective numerical optimization (https://github.com/P-N-Suganthan/2024-CEC).

Can be used for CEC 2017, CEC 2022 and CEC 2024 benchmarks

The algorithm code is in "L-SRTDE.cpp" file.

# Compilation and usage

Compilation is simple using gcc/g++:

g++ -std=c++11 -O3 L-SRTDE.cpp -o L-SRTDE.exe

or depending on hardware

g++ -std=c++11 -O3 -march=corei7-avx L-SRTDE.cpp -o L-SRTDE.exe

Please note that the compilation requires support of C++11 standard.

This will create L-SRTDE executable, available for running.

Data will be written to "L-NTADE_F#_D#.txt", where F and DIM are the function number and problem dimention.
