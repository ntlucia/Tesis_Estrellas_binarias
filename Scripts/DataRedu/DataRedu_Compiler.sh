#!/bin/bash

echo "#!/bin/bash

echo -e \"SPECTRUM\tRESOLUTION\tSNR\tIRF\" > \`echo \${PWD##*/}\`_INFO.dat

if [ -d *\"Fits\" ]
then
    cd ./\${PWD##*/}\ Fits
    rm -f *vo.fits
    rm -f *small.fits
    
    python $(pwd | sed 's/ /\\ /g')/con_to_ascii_v2_iSpec.py
    
    mv ./*dat ../
    cd ..

    #python $(pwd | sed 's/ /\\ /g')/iSRedu_B.py
    
    DIR=\"\${PWD##*/}\"
    
    for f in ./\${PWD##*/}\ Fits/*R*.fits; do
        SPECTRUM=\"\$(echo \"\$f\" | sed 's| ||g' | sed \"s|./\${DIR}Fits/Sci_||g\" | sed 's|.sp_ech_main.fits||g')\"
	
        MRESOL=\"\$(grep -oP --text \"MRESOL.{25}\" ./\${PWD##*/}\ Fits/*\${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')\"
	
        SNR=\"\$(grep -oP --text \"MEAN_SNR.{15}\" ./\${PWD##*/}\ Fits/*\${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')\"
	
        IRF=\"\$(grep -oP --text \"IRF_CORR.{20}\" ./\${PWD##*/}\ Fits/*\${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')\"
	
        echo -e \"\$SPECTRUM\t\$MRESOL\t\$SNR\t\$IRF\" >> \`echo \${PWD##*/}\`_INFO.dat
    done
    
else
    mkdir ./\${PWD##*/}\ Fits
    mv ./*fits ./\${PWD##*/}\ Fits       
    
    cd ./\${PWD##*/}\ Fits
    rm -f *vo.fits
    rm -f *small.fits
    
    python $(pwd | sed 's/ /\\ /g')/con_to_ascii_v2_iSpec.py
    
    mv ./*dat ../
    cd ..

    #python $(pwd | sed 's/ /\\ /g')/iSRedu_B.py

    DIR=\"\${PWD##*/}\"    

    for f in ./\${PWD##*/}\ Fits/*R*.fits; do
        SPECTRUM=\"\$(echo \"\$f\" | sed 's| ||g' | sed \"s|./\${DIR}Fits/Sci_||g\" | sed 's|.sp_ech_main.fits||g')\"
	
        MRESOL=\"\$(grep -oP --text \"MRESOL.{25}\" ./\${PWD##*/}\ Fits/*\${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')\"
	
        SNR=\"\$(grep -oP --text \"MEAN_SNR.{15}\" ./\${PWD##*/}\ Fits/*\${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')\"
	
        IRF=\"\$(grep -oP --text \"IRF_CORR.{20}\" ./\${PWD##*/}\ Fits/*\${SPECTRUM}*.fits | sed 's| ||g' | sed 's|^.*=||g' | sed 's|=||g' | sed 's/://g')\"
	
        echo -e \"\$SPECTRUM\t\$MRESOL\t\$SNR\t\$IRF\" >> \`echo \${PWD##*/}\`_INFO.dat
    done
fi"> DataRedu_iSpec.sh

echo "#!/bin/bash
SDIR=\"$(pwd | sed 's/ /\\ /g')\"
echo \"Script directory defined as: \$SDIR\"

EDIR=\"\$(pwd | sed 's/ /\\\ /g')\"
echo \"Execution directory defined as: \$EDIR\"

eval cd \$EDIR/

for dir in *; do [ -d \"\$dir\" ] && eval cp \$SDIR/DataRedu_iSpec.sh \"\$dir\" ; done
echo \"DataRedu_iSpec copied\"

eval cd \$EDIR/

for dir in *; do [ -d \"\$dir\" ] && eval cd \$dir && ./DataRedu_iSpec.sh && cd .. ; done
echo \"DataRedu_iSpec executed\"
" > DataRedu.sh

chmod 775 DataRedu_iSpec.sh
chmod 775 DataRedu.sh
