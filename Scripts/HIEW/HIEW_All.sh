#!/bin/bash
SDIR="/c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/HIEW"
echo "Script directory defined as: $SDIR"

EDIR="$(pwd | sed 's/ /\\ /g')"
echo "Execution directory defined as: $EDIR"

eval cd $EDIR/

for dir in *; do [ -d "$dir" ] && eval cp $SDIR/HIEW.sh "$dir" ; done
echo "HIEW  copied"

eval cd $EDIR/

for dir in *; do [ -d "$dir" ] && eval cd $dir && ./HIEW.sh && cd .. ; done
echo "HIEW executed"

