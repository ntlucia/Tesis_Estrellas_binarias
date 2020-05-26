# -------------------------------------------------------------------------------------------------------

NAME = system("echo ${PWD##*/}")
load "HIEW_Results.txt"

# -------------------------------------------------------------------------------------------------------

set xzeroaxis
set yzeroaxis
set xrange [3933.66-2:3933.66+2]
set offsets 0, 0, 0.006, 0
set xtics 0.5
set mxtics 5
set ytics 0.01
set mytics 2

set encoding iso_8859_1
set terminal pdf enhanced color font "DejaVuSans, 15" size 6,4.5
set output "HIEW_Plot_".NAME.".pdf"

set xlabel "Wavelength ({\305})"
set ylabel "Relative Intensity"

# -------------------------------------------------------------------------------------------------------

set label 1 "" at L_IMin_1,IMin_1 point pointtype 3 pointsize 2
set label 2 "" at L_IMax_1,IMax_1 point pointtype 3 pointsize 2
set label 3 "" at L_IMin_2,IMin_2 point pointtype 3 pointsize 2
set label 4 "" at L_IMax_2,IMax_2 point pointtype 3 pointsize 2
set label 5 "W_0=".system(sprintf("echo %0.2f {\261} %0.2f {\305}", FWHM, E_FWHM)) \
at graph(0.70,0.70), graph(0.90,0.90)
set label 6 system("echo ${PWD##*/}") at graph(0.075,0.075), graph(0.90,0.90)

# -------------------------------------------------------------------------------------------------------

set arrow from FWHM_LL,graph(0,0) to FWHM_LL,graph(1,1) nohead lt 2 lc 0 lw 1 dt 2
set arrow from FWHM_UL,graph(0,0) to FWHM_UL,graph(1,1) nohead lt 2 lc 0 lw 1 dt 2

set style rect fc lt -1 fs solid 0.25 noborder behind

set object rect from (FWHM_LL-ELFWHM_LL),graph(0,0) to (FWHM_LL+EUFWHM_LL),graph(1,1)
set object rect from (FWHM_UL-ELFWHM_UL),graph(0,0) to (FWHM_UL+EUFWHM_UL),graph(1,1)

# -------------------------------------------------------------------------------------------------------

plot Data_File using ($1*10):2 with histeps lc 0 lw 1 notitle, \
Spline_Out_File using 1:2 with lines lc 3 lw 2 notitle

#plot Data_File using ($1*10):2 with lines lc 0 lw 1 notitle, \
#Spline_Out_File using 1:2 with lines lc 3 lw 2 notitle
