#!/bin/bash

mode=$1
instance=$2
tolerance=$3

if [[ "$instance" =~ ^n.*$ ]]; then
    problem="smwtp"
else
    problem="arp"
fi



python3 fourier.py --mode $mode --problem $problem --instance instances/$problem/$instance --tolerance $tolerance > /dev/null

# python3 fourier.py --mode $mode --problem $problem --instance instances/$problem/$instance --tolerance ftolerance 2> $instance"_ftolerance.json"

# if [[ "$mode" != "YOR" ]]; then
#     python3 fourier.py --mode $mode --problem $problem --instance instances/$problem/$instance --tolerance r 2> $instance"_r.json"
# fi