#!/bin/bash

# This script must be run in the root directory of the github
test=$PWD/Scripts/06_Two_Thousand_Features_Dup/evaluate_imputation

# Clear output from previous runs
#rm ${test}/Error/*
#rm ${test}/Output/*
#rm ${test}/Results/*/*

# Run the jobs
for i in butterfly; do # Data Set
  for j in {1..3}; do # Imputation Method
    for k in 0; do # Seed
      echo "Rscript ${test}/job.R $i $j $k" | qsub -N "${i}_${j}_${k}" -l nodes=1:ppn=12,gres=localhd:1,vmem=90G,mem=90G,walltime=7:00:00:00 -o ${test}/Output -e ${test}/Error
      sleep 0.1
    done
    sleep 0.5
  done
done
