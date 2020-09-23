set term x11
set terminal postscript color enhanced font "Helvetica,20"
set output "Cat-Acc.ps"

set size square
set xlabel 'rho/10^{13} gcc'
#set ylabel 'A'
#set pointsize 1
set xrange [0:5]
#set yrange [29.5:35.5]
set mxtics 5
set mytics 5
set key samplen 1
plot \
 'Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat'    u ($1/1e13):4 w l lc 1 lt -1 lw 2 title 'A cell',\
 'Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat'    u ($1/1e13):5 w l lc 1 lt 4 lw 2 title 'A ion',\
 'Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat'    u ($1/1e13):6 w l lc 1 lt 2 lw 2 title 'Z',\
 'Crust_EOS_Cat_HZD-NV.dat'    u ($1/1e13):4 w l lc 0 lt -1 lw 2 title '',\
 'Crust_EOS_Cat_HZD-NV.dat'    u ($1/1e13):5 w l lc 0 lt 4 lw 2 title '',\
 'Crust_EOS_Cat_HZD-NV.dat'    u ($1/1e13):6 w l lc 0 lt 2 lw 2 title ''



