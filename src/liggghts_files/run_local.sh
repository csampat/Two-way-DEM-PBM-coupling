#!/bin/bash
bash rem.sh
mpiexec -n 4 mpi_liggghts -in in.final_2mm > output_2mm.txt
