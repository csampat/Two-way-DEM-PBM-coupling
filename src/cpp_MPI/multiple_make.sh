#!/bin/bash

for i in 64 128 256;
	do for j in 200;
		do vi "+:5s/.*/CORES="$i "+:wq" Makefile;
		   vi "+:6s/.*/DIAM="$j "+:wq" Makefile;
		   echo $i_$j;
		   # bash run.sh;
		   make clean;
		   make;
		   rm -f .Makefile.swp;
	   done
   done
