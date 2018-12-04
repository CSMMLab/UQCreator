#!/bin/bash
clear

BIN="../code/bin/UQCreator"

MIN_NCELLS=50
MAX_NCELLS=500
MIN_VAR=1
MAX_VAR=5
MIN_MOMENTS=3
MAX_MOMENTS=10
MIN_QUADPTS=10
MAX_QUADPTS=30

if [[ ! -d "BB" && ! -d "SG" && ! -d "SC" ]] ; then
    ./generateSettings.sh $MIN_NCELLS $MAX_NCELLS $MIN_VAR $MAX_VAR $MIN_MOMENTS $MAX_MOMENTS $MIN_QUADPTS $MAX_QUADPTS
fi

#BB_FILES=BB/input/settings*
#SG_FILES=SG/input/settings*
#SC_FILES=SC/input/settings*

#rm -f log

#for f in $(find $BB_FILES -iname "settings*" | sort -V ); do
#    echo -ne "Simulating $f...  \t"
#    ./$BIN -c $f >> log 2>&1
#    echo "DONE"
#done

#for f in $(find $SG_FILES -iname "settings*" | sort -V ); do
#    echo -ne "Simulating $f... \t"
#    ./$BIN -c $f >> log 2>&1
#    echo "DONE"
#done

#for f in $(find $SC_FILES -iname "settings*" | sort -V ); do
#    echo -ne "Simulating $f...   \t"
#    ./$BIN -c $f >> log 2>&1
#    echo "DONE"
#done

#echo -e "\n -> Finished!"
