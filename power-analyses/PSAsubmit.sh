#!/bin/bash
# File: PSAsubmit.sh

# Generate the sequences for $1 and $2 using seq command in bash
N=$(seq 50 25 250)
N_c=$(seq 30 10 50)
N_rep=50
iteration="A B C D E F"


# Iterate over the values for $1, $2 and iteration
for val1 in $N
do
    for val2 in $N_c
    do
        for letter in $iteration
        do
            # Submit the job with current values
            qsub PSAsim_wrapper.sh $val1 $val2 $N_rep $letter
        done
    done
done

