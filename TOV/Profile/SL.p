#output
set term x11
set terminal postscript enhanced color
set output "SL.ps"

##title
set title "EOS : Sly4/NV/HZD"


##plot n°1

##axis
set ylabel "Density g cm^{-3}"
set xlabel "Radius (km)"

##additional stuff

#set logscale y
#set xrange [1:1e14]
#set yrange [10.8:11.6]


plot \
'Prof_SL_Cat+NV_1.4_smooth.dat' using ($2/1e3):4 w l lc 1 lt 1 title '1.4 M_{sun}',\
'Prof_SL_Cat+NV_1.5_smooth.dat' using ($2/1e3):4 w l lc 2 lt 1 title '1.5 M_{sun}',\
'Prof_SL_Cat+NV_1.6_smooth.dat' using ($2/1e3):4 w l lc 3 lt 1 title '1.6 M_{sun}',\
'Prof_SL_Cat+NV_1.7_smooth.dat' using ($2/1e3):4 w l lc 4 lt 1 title '1.7 M_{sun}',\
'Prof_SL_Cat+NV_1.8_smooth.dat' using ($2/1e3):4 w l lc 5 lt 1 title '1.8 M_{sun}',\
'Prof_SL_Cat+NV_1.9_smooth.dat' using ($2/1e3):4 w l lc 6 lt 1 title '1.9 M_{sun}',\
'Prof_SL_Cat+NV_2.0_smooth.dat' using ($2/1e3):4 w l lc 8 lt 1 title '2.0 M_{sun}',\
'DURCA_threshold.dat' 		using 1:2 w l lc 7 lt 2 lw 2 title 'DURCA threshold'

