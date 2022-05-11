#!/bin/bash
# To be used with ligand_matcher.py

# prints out unique atoms to each COMBS group number
awk '{print $6,$2}' out.txt | sort -n | uniq |
# places each atom name on the same line next to corresponding group number
awk '{if(a!=$1) {a=$1; printf "\n%s%s",$0,FS} else {a=$1;$1="";printf $0 }} END {printf "\n"}' |
# replaces characters with commas to create a readable CSV file
sed '/^$/d' | sed 's/ /,/' | sed 's/ //' > sorted.csv
