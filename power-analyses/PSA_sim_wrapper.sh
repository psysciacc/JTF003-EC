#!/bin/bash
# File: PSAsim_wrapper.sh
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -q all.q
#$ -l mem_free=10G
#$ -m be
#$ -M ml19179@essex.ac.uk
./PSA_pwr_script_cluster.R $1 $2 $3 $4 $5
