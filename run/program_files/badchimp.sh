#!/bin/bash

### check number of args
if [[ $# -lt 3 ]]; then
    echo "Usage: "
    echo "       badchimp.sh [run-dir] [bin-type] [geo-file]"
    echo
    exit 1
fi


### set variables from arguments
DIR=$( echo $(dirname "${BASH_SOURCE[0]}") | awk -F '/' '{print $2}')   # directory of this script
RUN="$1"
BIN="$DIR/bin/$2/lb_ejah"
GEO="$DIR/geo/$3"
INP="$RUN/input.dat"
OPT="$RUN/options.dat"
BAS="$DIR/basis"
DAT="$DIR/data"

echo $DIR

### check if out folder exists, and increment if true
c=1
OUT="$RUN/out"
while [ -d "$OUT" ]; do
    OUT="$RUN/out-$((c++))"
done

### currently need this if running on a mac
MCA=""
if [[ $OSTYPE == *"darwin"* ]]; then
    MCA="--mca io romio314"
fi

### extract mpi vector from input.dat
NP=$(grep -hre 'num_proc' "$INP" | awk '{print $2*$3*$4}')

echo "mpirun ${MCA} -n ${NP} .${BIN} -basis ${BAS} -data ${DAT} -out ${OUT} -geo ${GEO} -input ${INP} -options ${OPT}"
mpirun ${MCA} -n ${NP} ${BIN} -basis ${BAS} -data ${DAT} -out ${OUT} -geo ${GEO} -input ${INP} -options ${OPT}


