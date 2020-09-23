
#output
set term x11
set terminal postscript enhanced color
set output "compa.ps"

##title
set title "Comparison between the two SL EOS : M=1.6 M_{sol}"


##plot n°1

##axis
set xlabel "Density g cm^{-3}"
set ylabel "Radius (km)"

##additional stuff

set logscale x
set xrange [1:1e14]
set yrange [10.8:11.6]


plot 'Prof_SL_Cat_1.6.dat' using 4:($2/1e3) title 'DH + HZD' with lines\
     ,'Prof_SL_Cat+NV_1.6.dat' using 4:($2/1e3) title 'DH + NV + HZD' with lines


##plot n°2

unset xrange
unset yrange
unset logscale xy

##axis
set xlabel "Density g cm^{-3} "
set ylabel "Enclosed Mass (M_{sol})"


set logscale x
set xrange [1:1e14]
set yrange [1.585:1.605]

plot 'Prof_SL_Cat_1.6.dat' using 4:6 title 'DH + HZD' with lines\
     ,'Prof_SL_Cat+NV_1.6.dat' using 4:6 title 'DH + NV + HZD' with lines

##plot n°3

unset xrange
unset yrange
unset logscale xy

##axis
set xlabel "Density g cm^{-3} "
set ylabel "Pressure (dyn cm^{-2}"

set logscale x
set logscale y
set xrange [1:1e14]
set yrange [1:1e33]

plot 'Prof_SL_Cat_1.6.dat' using 4:5 title 'DH + HZD'w l\
     ,'Prof_SL_Cat+NV_1.6.dat' using 4:5 title 'DH + NV + HZD' w l


