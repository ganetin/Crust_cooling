#   Rho	   Press	    nbar	    Ye	    Ymu	    Yn	    Yp	    Yla	    Ysm	    Ys0	    Ysp	Yxim	  
#    1      2		     3              4        5       6       7       8       9       10      11  12     

set term x11
set terminal postscript eps color enhanced font "Arial, 10" fontfile '/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb'
set output "Durca_BM165S_e.eps"

kf(x)=(x*3.*3.14159**2.)**(1/3.)
f0(x)=0.

set size 1,1
set origin 0,0
set multiplot layout 3,2 

set xlabel 'n_b [fm^{-3}]'
set ylabel 'Triangles'

set title 'DURCAe np'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($4*$3)+kf($6*$3)-kf($7*$3)) w l ls 1 lc 4 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($6*$3)+kf($7*$3)-kf($4*$3)) w l ls 2 lc 4 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($7*$3)+kf($4*$3)-kf($6*$3)) w l ls 3 lc 4 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCAe {/Symbol L}^0p'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($4*$3)+kf($8*$3)-kf($7*$3)) w l ls 1 lc 3 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($8*$3)+kf($7*$3)-kf($4*$3)) w l ls 2 lc 3 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($7*$3)+kf($4*$3)-kf($8*$3)) w l ls 3 lc 3 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCAe {/Symbol S}^-{/Symbol L}^0'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($4*$3)+kf($8*$3)-kf($9*$3)) w l ls 1 lc 8 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($8*$3)+kf($9*$3)-kf($4*$3)) w l ls 2 lc 8 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($9*$3)+kf($4*$3)-kf($8*$3)) w l ls 3 lc 8 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCAe {/Symbol X}^-{/Symbol L}^0'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($4*$3)+kf($12*$3)-kf($8*$3)) w l ls 1 lc 1 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($12*$3)+kf($8*$3)-kf($4*$3)) w l ls 2 lc 1 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($8*$3)+kf($4*$3)-kf($12*$3)) w l ls 3 lc 1 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCAe {/Symbol S}^-n'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($4*$3)+kf($6*$3)-kf($9*$3)) w l ls 1 lc 7 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($6*$3)+kf($9*$3)-kf($4*$3)) w l ls 2 lc 7 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($9*$3)+kf($4*$3)-kf($6*$3)) w l ls 3 lc 7 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''
unset multiplot


set term x11
set terminal postscript eps color enhanced font "Arial, 10" fontfile '/usr/share/texmf-texlive/fonts/type1/public/amsfonts/cm/cmsy10.pfb'
set output "Durca_BM165S_mu.eps"

kf(x)=(x*3.*3.14159**2.)**(1/3.)
f0(x)=0.

set size 1,1
set origin 0,0
set multiplot layout 3,2 

set xlabel 'n_b [fm^{-3}]'
set ylabel 'Triangles'

set title 'DURCA{/Symbol m} np'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($5*$3)+kf($6*$3)-kf($7*$3)) w l ls 1 lc 4 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($6*$3)+kf($7*$3)-kf($5*$3)) w l ls 2 lc 4 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($7*$3)+kf($5*$3)-kf($6*$3)) w l ls 3 lc 4 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCA{/Symbol m} {/Symbol L}^0p'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($5*$3)+kf($8*$3)-kf($7*$3)) w l ls 1 lc 3 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($8*$3)+kf($7*$3)-kf($5*$3)) w l ls 2 lc 3 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($7*$3)+kf($5*$3)-kf($8*$3)) w l ls 3 lc 3 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCA{/Symbol m} {/Symbol S}^-{/Symbol L}^0'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($5*$3)+kf($8*$3)-kf($9*$3)) w l ls 1 lc 8 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($8*$3)+kf($9*$3)-kf($5*$3)) w l ls 2 lc 8 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($9*$3)+kf($5*$3)-kf($8*$3)) w l ls 3 lc 8 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCA{/Symbol m} {/Symbol X}^-{/Symbol L}^0'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($5*$3)+kf($12*$3)-kf($8*$3)) w l ls 1 lc 1 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($12*$3)+kf($8*$3)-kf($5*$3)) w l ls 2 lc 1 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($8*$3)+kf($5*$3)-kf($12*$3)) w l ls 3 lc 1 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

set title 'DURCA{/Symbol m} {/Symbol S}^-n'
plot './BM165S_EOS_Cat_new_core.dat' u 3:(kf($5*$3)+kf($6*$3)-kf($9*$3)) w l ls 1 lc 7 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($6*$3)+kf($9*$3)-kf($5*$3)) w l ls 2 lc 7 title '',\
     './BM165S_EOS_Cat_new_core.dat' u 3:(kf($9*$3)+kf($5*$3)-kf($6*$3)) w l ls 3 lc 7 title '',\
     f0(x) w l ls 1 lc 7 lw 2 title ''

unset multiplot
