#!/bin/bash

# Remove all image files in the current directory
rm -f *.png 

# Loop from 0 to 10 inclusive
for i in {0..10}; do
    echo "Running for RHYTHM=$i"
    make gen RHYTHM=$i
    make run
done