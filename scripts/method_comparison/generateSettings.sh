#!/bin/bash

if [ "$#" -ne 8 ]; then
    echo "8 input parameters required!"
    exit
fi

#available methods
CLOSURES=(BB SG SC)
TIMESTEPPINGS=(EE SSPMS)

#store inputs
MIN_NCELLS=$1
MAX_NCELLS=$2
MIN_VAR=$3
MAX_VAR=$4
MIN_MOMENTS=$5
MAX_MOMENTS=$6
MIN_QUADPTS=$7
MAX_QUADPTS=$8

#replaces 1 to 2 in file 3
function replace {
    sed -i -e "s/${1}/${2}/g" ${3}
}

#create folders + copy and change files
for c in "${CLOSURES[@]}" ; do
    rm -rf $c/input
    mkdir -p $c $c/input $c/output
    for t in "${TIMESTEPPINGS[@]}" ; do
        for ((nc = MIN_NCELLS; nc <= MAX_NCELLS; nc=nc+50)); do
            CELLS=$(echo "$nc - 4" | bc)
            for ((v = MIN_VAR; v <= MAX_VAR; v++)); do
                for ((m = MIN_MOMENTS; m <= MAX_MOMENTS; m++)); do
                    for ((q = MIN_QUADPTS; q <= MAX_QUADPTS; q=q+5)); do
                        FILENAME="$c/input/settings_"
                        FILENAME+=$t
                        FILENAME+="_NC"
                        FILENAME+=$nc
                        FILENAME+="_VAR0"
                        FILENAME+=$v
                        FILENAME+="_MOM"
                        FILENAME+=$m
                        FILENAME+="_QUAD"
                        FILENAME+=$q
                        FILENAME+=".toml"
                        cp template.toml $FILENAME
                        replace CLOSURE $c $FILENAME
                        replace TIMESTEPPING $t $FILENAME
                        if [ $t == "EE"  ] ; then
                            if [ $c == "BB" ] ; then
                                replace CFLNUM 0.5 $FILENAME
                            else
                                replace CFLNUM 1.0 $FILENAME
                            fi
                        elif [ $t == "SSPMS" ] ; then
                            replace CFLNUM 0.33 $FILENAME
                        fi
                        replace NCELLS $(bc <<< "$nc-4") $FILENAME
                        replace VAR $(bc <<< "scale=1;$v/10") $FILENAME
                        replace MOMENTS $m $FILENAME
                        if [ $c == "SC" ] ; then
                            replace SC SG $FILENAME
                            replace QUADPTS $m $FILENAME
                            if [ $m -ne $q ] ; then
                                NEWFILENAME=$(echo $FILENAME | rev | cut -d'_' -f2- | rev)
                                NEWFILENAME+="_QUAD"
                                NEWFILENAME+=$m
                                NEWFILENAME+=".toml"
                                mv $FILENAME $NEWFILENAME
                            fi
                            break
                        else
                            replace QUADPTS $q $FILENAME
                        fi
                    done
                done
            done
        done
    done
done
