#!/bin/bash
SDIR="/c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/DataRedu"
echo "Script directory defined as: $SDIR"

EDIR="$(pwd | sed 's/ /\\ /g')"
echo "Execution directory defined as: $EDIR"

eval cd $EDIR/

for dir in *; do [ -d "$dir" ] && eval cp $SDIR/DataRedu_iSpec.sh "$dir" ; done
echo "DataRedu_iSpec copied"

eval cd $EDIR/

for dir in *; do [ -d "$dir" ] && eval cd $dir && ./DataRedu_iSpec.sh && cd .. ; done
echo "DataRedu_iSpec executed"

