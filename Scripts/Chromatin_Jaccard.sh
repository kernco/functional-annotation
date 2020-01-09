#!/bin/bash

echo -n "State\t"
for pre1 in $2; do
    for pre2 in $2; do
        echo -n "${pre1}_${pre2}\t"
    done
done
echo

for state in $(seq 1 $1); do
    echo -n $state "\t"
    for pre1 in $2; do
        for pre2 in $2; do
            rm temp1 temp2
            grep E$state$ ${pre1}_$3 | grep "^1\s" > temp1
            grep E$state$ ${pre2}_$3 | grep "^1\s" > temp2
            echo -n `bedtools jaccard -a temp1 -b temp2 | tail -n1 | cut -f 3` "\t"
        done
    done
    echo
done
