set term x11
set terminal postscript color enhanced font "Helvetica,20"
set output "Meff_BSK21.ps"

set size square
set xlabel 'rho/10^{15} gcc'
#set ylabel 'meff/m'
#set pointsize 1
#set xrange [-13:-8]
#set yrange [29.5:35.5]
set mxtics 5
set mytics 5
set key samplen 2
plot \
 'BSK21_EOS_Acc_Fe.dat'    u ($1/1e15):12 w l lc 0 lt 1 lw 2 title 'm_{p}^{eff}',\
 'BSK21_EOS_Acc_Fe.dat'    u ($1/1e15):13 w l lc 1 lt 1 lw 2 title 'm_{n}^{eff}'
