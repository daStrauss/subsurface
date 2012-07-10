#!/bin/bash

for ((i=172;i<=200;++i));
do
    cd "trial$i"
    rm *.nts
    cd ..
done