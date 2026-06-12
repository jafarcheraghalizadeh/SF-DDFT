#!/bin/bash

# Compile the Fortran program
gfortran -O3 extract_last_sample.f90 -o extract_last_sample


# Check if the compilation was successful
if [ $? -ne 0 ]; then
    echo "Fortran compilation failed!"
    exit 1
fi

# Loop over files sample1.dat to sample20.dat
for i in {1..20}; do
    # Set input file name
    
    input_file="sample${i}.dat"
    
    # Set output file name for the last sample
    output_file="last_sample${i}.dat"
    
    # Run the Fortran program (compiled as `extract_last_sample`)
    nohup ./extract_last_sample "$input_file" "$output_file" &
    
    # Check if the program ran successfully
    if [ $? -eq 0 ]; then
        echo "Processed $input_file and wrote last sample to $output_file"
    else
        echo "Error processing $input_file"
    fi
done

