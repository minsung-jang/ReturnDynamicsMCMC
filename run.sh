#!/bin/bash
cmake -E make_directory build
cd build
cmake ..
make
./ReturnDynamicsMCMC
cd ../
#cd ../plot_script
#python3 plot_script.py
