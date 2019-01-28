#!/bin/sh
#for i in {1..}
declare -a arr=("coi1" "npr1")
for j in "${arr[@]}"
do
  for i in {1..23956}
  #in reality, columns 6:23961
  do
    echo "genotype $j"
    echo "Looping ... number $i"
  #choosing to run with only kmat 2
    ~/gemma/bin/gemma -bfile $j/02_GEMMA/binMAF20NA10 -k $j/03_kmat/binMAF20NA10_kmat2_pheno1.sXX.txt -n $i -miss 0.1 -maf 0.2 -lm 4 -o "${j}"_MAF20NA10"_${i}"
  done
done
