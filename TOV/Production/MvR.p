
#output
set term x11
set terminal postscript enhanced color
set output "MvR.ps"

##title
set title "Comparison between the two SL EOS"



##axis
set xlabel "Radius (km)"
set ylabel "Enclosed Mass (M_{sol})"

##additional stuff

#set logscale x
set xrange [7:25]
#set yrange [10.8:11.6]


plot 'Prod_SL_EOS_Cat.dat' using 2:3 title 'DH + HZD' with lines\
     ,'Prod_SL_EOS_Cat+NV.dat' using 2:3 title 'DH + NV + HZD' with lines

