#!/bin/bash

# Compile the Fortran code (if not already compiled)
make clean
make main_blend

# Loop to generate files a1.dat, a2.dat, etc.
for i in {1..20}; do
    filename="sample${i}.dat"
    filename_dens="last_sample${i}.dat"
    outputfile="run${i}.dat"
    echo "Generating file: $filename"
    echo "PID:" $!
    # Run the Fortran program with the filename as an argument
    nohup ./main_blend.exe $filename $filename_dens > $outputfile 2>&1 & 
   #./main_multiblock.exe $filename $filename_dens
     sleep 2
done

