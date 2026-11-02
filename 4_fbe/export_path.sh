#!/bin/bash
# Script to set LD_LIBRARY_PATH

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/path/to/lib/:/scratch/path/to/lib/lib/

echo "LD_LIBRARY_PATH set to: $LD_LIBRARY_PATH"

export PATH=/scratch/path/to/lib/openmpi/bin:$PATH; export LD_LIBRARY_PATH=/path/to/lib/openmpi/lib:$LD_LIBRARY_PATH

