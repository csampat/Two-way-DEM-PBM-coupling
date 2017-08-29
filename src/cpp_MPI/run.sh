#!/bin/bash

rm -f *.txt

make clean

make -j 7

echo "compilation & linkin done"

echo "Running model.out"

#export OMP_NUM_THREADS=2

mpirun -n 2 ./model.out 128 200 > output_test.txt

echo "DONE"
