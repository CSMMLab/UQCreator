#!/bin/bash

AOAS=(0.0 0.25 0.5 0.75 1.0 1.25 1.5 1.75 2.0 2.25 2.5)

function replace {
    sed -i -e "s/${1}/${2}/g" ${3}
}

CORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')

rm -rf results
mkdir results

for ((i = 0; i<${#AOAS[@]}; i++)); do
    AOA=${AOAS[$i]}
    rm -rf $AOA
    mkdir -p $AOA
    cp template $AOA/$AOA.cfg
    replace AOA_VAL $AOA $AOA/$AOA.cfg
    replace OUTPUTFILENAME AOA_$AOA $AOA/$AOA.cfg
done

for f in $(find -iname "*.cfg" | sort -g); do
    FILE=${f##*/}
    DIR=${f%/*}
    echo -ne "Simulating AOA '$FILE'...   \t"
    cd $DIR
    mpirun -n $CORES SU2_CFD $FILE >> stdout 2>&1
    SU2_SOL $FILE >> stdout 2>&1
    mv AOA* ../results
    cd ..
    rm -rf $DIR
    echo "[DONE]"
done
