#!/bin/bash
g++ -o HIEW.out HIEW.cpp

echo "#!/bin/bash

if [ ! -f HIEW_Config.cfg ]
then
echo \"FILE=A_\${PWD##*/}_B.dat
OUTR=HIEW_Results.txt
OUTD=HIEW_Data.txt
OUTS=HIEW_Spline.txt
LMIN=3931.66
LMAX=3935.66
LLIN=3933.66
RANG=0.2
NSPL=200\" > HIEW_Config.cfg
fi

$(pwd | sed 's/ /\\ /g')/HIEW.out

gnuplot -p -e \"load '\"$(pwd | sed 's/ /\\ /g')/Plot_Script_HIEW.plt\"'\"

if [ -d *\"HIEW Results\" ]
then
    cd ./\${PWD##*/}\ HIEW\ Results
    mv ../HIEW_*.txt ./
    mv ../HIEW_*.pdf ./
    cd ..
else
    mkdir ./\${PWD##*/}\ HIEW\ Results
    cd ./\${PWD##*/}\ HIEW\ Results
    mv ../HIEW_*.txt ./
    mv ../HIEW_*.pdf ./
    cd ..
fi
" > HIEW.sh

echo "#!/bin/bash
SDIR=\"$(pwd | sed 's/ /\\ /g')\"
echo \"Script directory defined as: \$SDIR\"

EDIR=\"\$(pwd | sed 's/ /\\\ /g')\"
echo \"Execution directory defined as: \$EDIR\"

eval cd \$EDIR/

for dir in *; do [ -d \"\$dir\" ] && eval cp \$SDIR/HIEW.sh \"\$dir\" ; done
echo \"HIEW  copied\"

eval cd \$EDIR/

for dir in *; do [ -d \"\$dir\" ] && eval cd \$dir && ./HIEW.sh && cd .. ; done
echo \"HIEW executed\"
" > HIEW_All.sh

chmod 775 HIEW.sh
chmod 775 HIEW_All.sh
