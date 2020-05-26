#!/bin/bash

echo -e "SPECTRUM\tRESOLUTION\tSNR\tIRF" > `echo ${PWD##*/}`_INFO.dat

if [ -d *"Fits" ]
then
    cd ./${PWD##*/}\ Fits
    rm -f *vo.fits
    rm -f *small.fits
    
    python /c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/DataRedu/con_to_ascii_v2_iSpec.py
    
    mv ./*dat ../
    cd ..

    #python /c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/DataRedu/iSRedu_B.py
    
    DIR="${PWD##*/}"
    
    for f in ./${PWD##*/}\ Fits/*R*.fits; do
        SPECTRUM="$(echo "$f" | sed 's| ||g' | sed "s|./${DIR}Fits/Sci_||g" | sed 's|.sp_ech_main.fits||g')"
	
        MRESOL="$(grep -oP --text "MRESOL.{25}" ./${PWD##*/}\ Fits/*${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')"
	
        SNR="$(grep -oP --text "MEAN_SNR.{15}" ./${PWD##*/}\ Fits/*${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')"
	
        IRF="$(grep -oP --text "IRF_CORR.{20}" ./${PWD##*/}\ Fits/*${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')"
	
        echo -e "$SPECTRUM\t$MRESOL\t$SNR\t$IRF" >> `echo ${PWD##*/}`_INFO.dat
    done
    
else
    mkdir ./${PWD##*/}\ Fits
    mv ./*fits ./${PWD##*/}\ Fits       
    
    cd ./${PWD##*/}\ Fits
    rm -f *vo.fits
    rm -f *small.fits
    
    python /c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/DataRedu/con_to_ascii_v2_iSpec.py
    
    mv ./*dat ../
    cd ..

    #python /c/Users/lucia/Documents/Semestre\ 8\ y9/Trabajo\ grado/Scripts/DataRedu/iSRedu_B.py

    DIR="${PWD##*/}"    

    for f in ./${PWD##*/}\ Fits/*R*.fits; do
        SPECTRUM="$(echo "$f" | sed 's| ||g' | sed "s|./${DIR}Fits/Sci_||g" | sed 's|.sp_ech_main.fits||g')"
	
        MRESOL="$(grep -oP --text "MRESOL.{25}" ./${PWD##*/}\ Fits/*${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')"
	
        SNR="$(grep -oP --text "MEAN_SNR.{15}" ./${PWD##*/}\ Fits/*${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')"
	
        IRF="$(grep -oP --text "IRF_CORR.{20}" ./${PWD##*/}\ Fits/*${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')"
	
        echo -e "$SPECTRUM\t$MRESOL\t$SNR\t$IRF" >> `echo ${PWD##*/}`_INFO.dat
    done
fi
