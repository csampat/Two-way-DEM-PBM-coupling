#!/bin/bash
#SBATCH -J LIGGGHTSS_o3_100M
#SBATCH -o LIGGGHTSScalingTest.out
#SBATCH -e LIGGGHTSScalingTest.err
#SBATCH -n256 -N4
#SBATCH -p normal
#SBATCH -t 12:00:00
#SBATCH -A TG-MCB090174
#SBATCH --mail-user=crs234@scarletmail.rutgers.edu
#SBATCH --mail-type=end

#export OMP_NUM_THREADS=2 

#ulimit -c unlimited
ibrun -n 256 -o 0 tacc_affinity liggghts_micstam -in in.final_2mm >output_2mm.out

#ibrun tacc_affinity /work/04085/tg833843/liggghts_work/liggghts_builds_for_franklin_stamp/LIGGGHTS/LIGGGHTS-PUBLIC/src/./lmp_stampede_vtk -in in.cluster_5000000
