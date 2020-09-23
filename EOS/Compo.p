set term x11
set terminal postscript color enhanced font "Helvetica,20"# fontfile '/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb'

set output "compo.eps"
set colors classic
set xrange [0.2:1.2]
set yrange [1e-3:1]
set xtics 0.5
set mxtics 5  
set logscale y
set size square
set key at 1.05,4e-2

set xlabel "n_{B} [fm^{-3}]"
set ylabel "Y_i"
plot \
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($4>=0?$4:1/0)   w l lc 1  lt 1 linewidth 2  title 'e',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($5>=0?$5:1/0)   w l lc 2  lt 1 linewidth 2  title '{/Symbol m}',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($6>=0?$6:1/0)   w l lc 3  lt 1 linewidth 2  title 'n',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($7>=0?$7:1/0)   w l lc 4  lt 1 linewidth 2  title 'p',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($8>=0?$8:1/0)   w l lc 9  lt 1 linewidth 2  title '{/Symbol L}',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($9>=0?$9:1/0)   w l lc 1  lt 2 linewidth 2  title '{/Symbol S}^-',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($10>=0?$10:1/0) w l lc 2  lt 2 linewidth 2  title '{/Symbol S}^0',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($11>=0?$11:1/0) w l lc 3  lt 2 linewidth 2  title '{/Symbol S}^+',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($12>=0?$12:1/0) w l lc 4  lt 2 linewidth 2  title '{/Symbol X}^-',\
 'NH1_EOS_Acc_Fe_new_2.dat'  u 3:($13>=0?$13:1/0) w l lc 9  lt 2 linewidth 2  title '{/Symbol X}^0',\
 'G85_Yold_Acc_Fe.dat'  u 3:($4>=0?$4:1/0)   w l lc 1  lt 3 linewidth 2  title 'e',\
 'G85_Yold_Acc_Fe.dat'  u 3:($5>=0?$5:1/0)   w l lc 2  lt 3 linewidth 2  title '{/Symbol m}',\
 'G85_Yold_Acc_Fe.dat'  u 3:($6>=0?$6:1/0)   w l lc 3  lt 3 linewidth 2  title 'n',\
 'G85_Yold_Acc_Fe.dat'  u 3:($7>=0?$7:1/0)   w l lc 4  lt 3 linewidth 2  title 'p',\
 'G85_Yold_Acc_Fe.dat'  u 3:($8>=0?$8:1/0)   w l lc 9  lt 3 linewidth 2  title '{/Symbol L}',\
 'G85_Yold_Acc_Fe.dat'  u 3:($9>=0?$9:1/0)   w l lc 1  lt 4 linewidth 2  title '{/Symbol S}^-',\
 'G85_Yold_Acc_Fe.dat'  u 3:($10>=0?$10:1/0) w l lc 2  lt 4 linewidth 2  title '{/Symbol S}^0',\
 'G85_Yold_Acc_Fe.dat'  u 3:($11>=0?$11:1/0) w l lc 3  lt 4 linewidth 2  title '{/Symbol S}^+',\
 'G85_Yold_Acc_Fe.dat'  u 3:($12>=0?$12:1/0) w l lc 4  lt 4 linewidth 2  title '{/Symbol X}^-',\
 'G85_Yold_Acc_Fe.dat'  u 3:($13>=0?$13:1/0) w l lc 9  lt 4 linewidth 2  title '{/Symbol X}^0'#,\

# 'BM.dat' u 1:2 w l lc 7 lt 1 lw 3 title ''
#  Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	Yxi0
#   1        2               3               4       5       6      7        8       9       10      11  12      13


