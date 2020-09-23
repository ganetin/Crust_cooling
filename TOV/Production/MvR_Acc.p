
#output
set term x11
set terminal postscript enhanced color
set output "MvR.eps"

##title



##axis
set xlabel "Radius (km)"
set ylabel "Mass (M_{sol})"

##additional stuff

#set logscale x
set xrange [7:25]
#set yrange [10.8:11.6]


plot 'prod_APR_EOS_Acc_Fe.dat' using 2:3 title '' lc 7 with lines

