#!/bin/bash

make clean
make vcg

flag_key=1
coll_key=0

# Directory for output data
output_dir="output_data"
input_dir="../data_vis"

# Make sure the output directory exists
mkdir -p $output_dir

./vcg $input_dir/1000004080.fits inputparameters $output_dir/GV_data_1 $flag_key $coll_key 1 $output_dir/sampling

#for realization_number in {2..20}
#do
#    ./vcg $input_dir/mwa_${realization_number}.fits inputparameters $output_dir/GV_data_${realization_number} $flag_key $coll_key 0 ph
#    echo "############${realization_number} completed.###########################"
#done

echo "VCG Done."

##./vcg ../2b_visfitsgrid/data/Usingle_chan_${realization_number}.fits inputparameters $output_dir/UGV_data_${realization_number}

