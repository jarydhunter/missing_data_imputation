#!/bin/bash

home=$PWD
test=${home}Scripts/06_Two_Thousand_Features_Dup/cluster_complete_data

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in 4; do # Data Set
  echo "${test}/job.R $i" | qsub -N "${i}" -l nodes=1:ppn=35,gres=localhd:1,vmem=70G,mem=70G,walltime=4:00:00:00 -o ${test}/Output -e ${test}/Error
  sleep 0.1
done
