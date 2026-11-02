#!/bin/bash


make amber_v2

./amber_v2 1
echo "amber_v2 done"

#for realization_number in {101..250}
#do
#    ./mabs_final_NF ${realization_number}
#    echo "mabs_NF ${realization_number} done"
#done



#echo "All realizations completed."


