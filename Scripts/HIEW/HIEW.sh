#!/bin/bash

if [ ! -f HIEW_Config.cfg ]
then
echo "FILE=A_${PWD##*/}_B.dat
OUTR=HIEW_Results.txt
OUTD=HIEW_Data.txt
OUTS=HIEW_Spline.txt
LMIN=3931.66
LMAX=3935.66
LLIN=3933.66
RANG=0.2
NSPL=200" > HIEW_Config.cfg
fi

/c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/HIEW/HIEW.out

gnuplot -p -e "load '"/c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/HIEW/Plot_Script_HIEW.plt"'"

if [ -d *"HIEW Results" ]
then
    cd ./${PWD##*/}\ HIEW\ Results
    mv ../HIEW_*.txt ./
    mv ../HIEW_*.pdf ./
    cd ..
else
    mkdir ./${PWD##*/}\ HIEW\ Results
    cd ./${PWD##*/}\ HIEW\ Results
    mv ../HIEW_*.txt ./
    mv ../HIEW_*.pdf ./
    cd ..
fi

