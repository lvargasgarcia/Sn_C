#!/bin/bash

PROBLEM=$1

FILE=instances_$PROBLEM.txt

for i in {5,6,7}; do
    mkdir -p ./results/YKR/$PROBLEM/$i
    mkdir -p ./results/YOR/$PROBLEM/$i
    mkdir -p ./results/YSR/$PROBLEM/$i
done

while read line; do
    python3 fourier.py --mode YKR --problem $PROBLEM --instance $line
    python3 fourier.py --mode YOR --problem $PROBLEM --instance $line
    python3 fourier.py --mode YSR --problem $PROBLEM --instance $line 
done < $FILE