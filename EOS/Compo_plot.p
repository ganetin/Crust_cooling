set term x11
set terminal postscript color enhanced font "Helvetica,20"# fontfile '/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb'

set output "BM165_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.5
set mxtics 5  
set logscale y
set size square
set key at 1.05,4e-2

set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($4>=0?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($5>=0?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($6>=0?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($7>=0?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($8>=0?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($9>=0?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($10>=0?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($11>=0?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($12>=0?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'BM165S_EOS_Acc_Fe_new_2.dat'  u 3:($13>=0?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'BM.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13


set output "GM1A_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.5
set mxtics 5  
set logscale y
set size square
unset key
set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'GM1A_EOS_Acc_Fe.dat'  u 3:($4>=5e-4?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($5>=5e-4?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($6>=5e-4?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($7>=5e-4?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($8>=5e-4?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($9>=5e-4?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($10>=5e-4?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($11>=5e-4?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($12>=5e-4?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($13>=5e-4?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'GM1A.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13

set output "GM1B_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.5
set mxtics 5  
set logscale y
set size square
unset key
set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'GM1B_EOS_Acc_Fe.dat'  u 3:($4>=5e-4?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($5>=5e-4?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($6>=5e-4?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($7>=5e-4?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($8>=5e-4?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($9>=5e-4?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($10>=5e-4?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($11>=5e-4?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($12>=5e-4?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($13>=5e-4?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'GM1B.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13

set output "TM1C_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.5
set mxtics 5  
set logscale y
set size square
unset key
set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'TM1C_EOS_Acc_Fe.dat'  u 3:($4>=5e-4?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($5>=5e-4?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($6>=5e-4?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($7>=5e-4?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($8>=5e-4?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($9>=5e-4?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($10>=5e-4?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($11>=5e-4?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($12>=5e-4?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($13>=5e-4?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'TM1C.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13


set output "z0.0_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.4
set mxtics 2  
set logscale y
set size square
unset key
set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'z0.0_EOS_Acc_Fe.dat'  u 3:($4>=0?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($5>=0?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($6>=0?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($7>=0?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($8>=0?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($9>=0?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($10>=0?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($11>=0?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($12>=0?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($13>=0?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'z0.0.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13


set output "z0.8_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.4
set mxtics 2  
set logscale y
set size square
unset key
set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'z0.8_EOS_Acc_Fe.dat'  u 3:($4>=0?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($5>=0?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($6>=0?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($7>=0?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($8>=0?$8:1/0)   w l lc 5  lt 1 linewidth 2  title '{/Symbol L}',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($9>=0?$9:1/0)   w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($10>=0?$10:1/0) w l lc 7  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($11>=0?$11:1/0) w l lc 8  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($12>=0?$12:1/0) w l lc 9  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($13>=0?$13:1/0) w l lc 10 lt 2 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13

set output "zsu6_compo.eps"

set xrange [0.1:1.2]
set yrange [1e-3:1]
set xtics 0.5
set mxtics 5  
set logscale y
set size square
unset key
set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'zsu6_EOS_Acc_Fe.dat'  u 3:($4>=0?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($5>=0?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($6>=0?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($7>=0?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($8>=0?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($9>=0?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($10>=0?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($11>=0?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($12>=0?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($13>=0?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'zsu6.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13


set output "GM1A_meff.eps"

set xrange [0.05:1.6]
set yrange [0.35:0.95]
set xtics 0.5
set mxtics 5  
set ytics 0.1
set mytics 2 

unset logscale y
set size square

set xlabel "n_{B} [fm^{-3}]"
set ylabel "m*_i/m_i"
plot \
 'GM1A_EOS_Acc_Fe.dat'  u 3:($14)  w l lc 1  lt 1 linewidth 2  title 'p',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($15)  w l lc 2  lt 1 linewidth 2  title 'n',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($8 >5e-4?$16:1/0)  w l lc 3  lt 1 linewidth 2  title '{/Symbol L}',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($9 >5e-4?$17:1/0)  w l lc 4  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($10>5e-4?$18:1/0)  w l lc 5  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($11>5e-4?$19:1/0)  w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($12>5e-4?$20:1/0)  w l lc 7  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'GM1A_EOS_Acc_Fe.dat'  u 3:($13>5e-4?$21:1/0)  w l lc 8  lt 1 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13         14      15
#16         17                18            19            20        21

set output "GM1B_meff.eps"

set xrange [0.05:1.2]
set yrange [0.35:0.95]
set xtics 0.5
set mxtics 5  
set ytics 0.1
set mytics 2 

unset logscale y
set size square

set xlabel "n_{B} [fm^{-3}]"
set ylabel "m*_i/m_i"
plot \
 'GM1B_EOS_Acc_Fe.dat'  u 3:($14)  w l lc 1  lt 1 linewidth 2  title 'p',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($15)  w l lc 2  lt 1 linewidth 2  title 'n',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($8 >5e-4?$16:1/0)  w l lc 3  lt 1 linewidth 2  title '{/Symbol L}',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($9 >5e-4?$17:1/0)  w l lc 4  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($10>5e-4?$18:1/0)  w l lc 5  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($11>5e-4?$19:1/0)  w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($12>5e-4?$20:1/0)  w l lc 7  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'GM1B_EOS_Acc_Fe.dat'  u 3:($13>5e-4?$21:1/0)  w l lc 8  lt 1 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13         14      15
#16         17                18            19            20        21

set output "TM1C_meff.eps"

set xrange [0.05:1.7]
set yrange [0.35:0.95]
set xtics 0.5
set mxtics 5  
set ytics 0.1
set mytics 2 

unset logscale y
set size square

set xlabel "n_{B} [fm^{-3}]"
set ylabel "m*_i/m_i"
plot \
 'TM1C_EOS_Acc_Fe.dat'  u 3:($14)  w l lc 1  lt 1 linewidth 2  title 'p',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($15)  w l lc 2  lt 1 linewidth 2  title 'n',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($8 >5e-4?$16:1/0)  w l lc 3  lt 1 linewidth 2  title '{/Symbol L}',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($9 >5e-4?$17:1/0)  w l lc 4  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($10>5e-4?$18:1/0)  w l lc 5  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($11>5e-4?$19:1/0)  w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($12>5e-4?$20:1/0)  w l lc 7  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'TM1C_EOS_Acc_Fe.dat'  u 3:($13>5e-4?$21:1/0)  w l lc 8  lt 1 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13         14      15
#16         17                18            19            20        21

set output "z0.0_meff.eps"

set xrange [0.0:1.2]
set yrange [0:1]
set xtics 0.4
set mxtics 2  
set ytics 0.1
set mytics 2 

unset logscale y
set size square

set xlabel "n_{B} [fm^{-3}]"
set ylabel "m*_i/m_i"
plot \
 'z0.0_EOS_Acc_Fe.dat'  u 3:($14)  w l lc 1  lt 1 linewidth 2  title 'p',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($15)  w l lc 2  lt 1 linewidth 2  title 'n',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($8 >5e-4?$16:1/0)  w l lc 3  lt 1 linewidth 2  title '{/Symbol L}',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($9 >5e-4?$17:1/0)  w l lc 4  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($10>5e-4?$18:1/0)  w l lc 5  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($11>5e-4?$19:1/0)  w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($12>5e-4?$20:1/0)  w l lc 7  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'z0.0_EOS_Acc_Fe.dat'  u 3:($13>5e-4?$21:1/0)  w l lc 8  lt 1 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13         14      15
#16         17                18            19            20        21


set output "z0.8_meff.eps"

set xrange [0.0:1.2]
set yrange [0:1]
set xtics 0.4
set mxtics 2  
set ytics 0.1
set mytics 2 

unset logscale y
set size square

set xlabel "n_{B} [fm^{-3}]"
set ylabel "m*_i/m_i"
plot \
 'z0.8_EOS_Acc_Fe.dat'  u 3:($14)  w l lc 1  lt 1 linewidth 2  title 'p',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($15)  w l lc 2  lt 1 linewidth 2  title 'n',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($8 >5e-4?$16:1/0)  w l lc 3  lt 1 linewidth 2  title '{/Symbol L}',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($9 >5e-4?$17:1/0)  w l lc 4  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($10>5e-4?$18:1/0)  w l lc 5  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($11>5e-4?$19:1/0)  w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($12>5e-4?$20:1/0)  w l lc 7  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'z0.8_EOS_Acc_Fe.dat'  u 3:($13>5e-4?$21:1/0)  w l lc 8  lt 1 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13         14      15
#16         17                18            19            20        21


set output "zsu6_meff.eps"

set xrange [0.0:1.2]
set yrange [0:1]
set xtics 0.4
set mxtics 2  
set ytics 0.1
set mytics 2 

unset logscale y
set size square

set xlabel "n_{B} [fm^{-3}]"
set ylabel "m*_i/m_i"
plot \
 'zsu6_EOS_Acc_Fe.dat'  u 3:($14)  w l lc 1  lt 1 linewidth 2  title 'p',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($15)  w l lc 2  lt 1 linewidth 2  title 'n',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($8 >5e-4?$16:1/0)  w l lc 3  lt 1 linewidth 2  title '{/Symbol L}',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($9 >5e-4?$17:1/0)  w l lc 4  lt 1 linewidth 2  title '{/Symbol S}^-',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($10>5e-4?$18:1/0)  w l lc 5  lt 1 linewidth 2  title '{/Symbol S}^0',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($11>5e-4?$19:1/0)  w l lc 6  lt 1 linewidth 2  title '{/Symbol S}^+',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($12>5e-4?$20:1/0)  w l lc 7  lt 1 linewidth 2  title '{/Symbol X}^-',\
 'zsu6_EOS_Acc_Fe.dat'  u 3:($13>5e-4?$21:1/0)  w l lc 8  lt 1 linewidth 2  title '{/Symbol X}^0'
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0	   mstp	   mstn	   mstla	   mstsm	   msts0	   mstsp	mstxim	mstxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13         14      15
#16         17                18            19            20        21
