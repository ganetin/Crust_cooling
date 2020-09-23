#set term x11
#set terminal postscript color enhanced 
#set output "EOS.ps"

unset logscale x
unset logscale y
#set logscale x
#set logscale y

#set xrange [11.8:13.6]
#set yrange [0:1300]
set title 'EOS in the crust'
set ylabel "A_{cell}"
set xlabel "log10 ({/Symbol r} [g cm^{-3}])"

plot 'Crust_EOS_Acc_Fe_pycno_HZ2008.dat' u (log10($1)):4  w lp title 'EOS old','Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat' u (log10($1)):4  w lp title 'EOS new'
 
#set xrange [11.8:13.6]
#set yrange [35:140]
set title 'EOS in the crust'
set ylabel "A_{ion}"
set xlabel "log10 ({/Symbol r} [g cm^{-3}])"

plot 'Crust_EOS_Acc_Fe_pycno_HZ2008.dat' u (log10($1)):5 w lp title 'EOS old','Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat' u (log10($1)):5 w lp title 'EOS new'

#set xrange [11.8:13.6]
#set yrange [10:32.5]
set title 'EOS in the crust'
set ylabel "Z_{ion}"
set xlabel "log10 ({/Symbol r} [g cm^{-3}])"

plot 'Crust_EOS_Acc_Fe_pycno_HZ2008.dat' u (log10($1)):6 w lp title 'EOS old','Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat' u (log10($1)):6 w lp title 'EOS new'

#set xrange [11.8:13.6]
#set yrange [0:1100]
set title 'EOS in the crust'
set ylabel "A_{cell}-A_{ion}"
set xlabel "log10 ({/Symbol r} [g cm^{-3}])"

plot 'Crust_EOS_Acc_Fe_pycno_HZ2008.dat' u (log10($1)):($4-$5) w lp title 'EOS old','Crust_EOS_Acc_Fe_pycno_HZ2008_3.dat' u (log10($1)):($4-$5) w lp title 'EOS new'


