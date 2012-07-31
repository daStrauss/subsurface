#!/bin/bash

for ((i=0;i<=199;++i));
do
    cd "trial$i"
    let b=$i+1
    for ((j=0;j<=15;++j));
    do
	mv "notes"$j"_"$b".nts" "notes"$j"_"$i".nts"
	mv "pout_contrastX"$b".mat" "pout_contrastX"$i".mat"
    done
    cd ..
done