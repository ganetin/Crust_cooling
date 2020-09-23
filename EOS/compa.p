set term x11
set terminal postscript eps color enhanced font "Helvetica,20"
set output "compo.eps"

#set logscale x
set logscale y
#unset xtics 
#unset ytics
#unset mxtics  
#unset mytics

set xrange [1.5:10]
set yrange [1e-5:100]

#set xtics 50
#set xtics 100
#set ytics 0.2
#set ytics 0.1
#set mxtics 10  
#set mytics 10
set mxtics 10  

# position of the legend
set key top right
set key nobox

#set size 0.65,0.65

set xlabel "{/Symbol r} [10^{14} g cm^{-3}]"
set ylabel "Y_i"

set format y "10^{%L}"


plot \
'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):4 w l lc 1 lt 1 lw 2 title 'Y_e',\
'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):5 w l lc 2 lt 1 lw 2 title "Y_{/Symbol m}",\
'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):6 w l lc 3 lt 1 lw 2 title 'Y_n',\
'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):7 w l lc 4 lt 1 lw 2 title 'Y_p',\
'APR_EOS_Cat.dat'            u ($1/1e14):4 w l lc 1 lt 2 lw 2 title '',\
'APR_EOS_Cat.dat'            u ($1/1e14):5 w l lc 2 lt 2 lw 2 title "",\
'APR_EOS_Cat.dat'   	     u ($1/1e14):6 w l lc 3 lt 2 lw 2 title '',\
'APR_EOS_Cat.dat'  	     u ($1/1e14):7 w l lc 4 lt 2 lw 2 title ''

#'APR_EOS_Cat.dat'            u ($1/1e14):4 w l lc 1 lt 1 lw 2 title 'APR - Ye',\
#'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):4 w l lc 1 lt 2 lw 2 title 'SL - Ye',\
#'APR_EOS_Cat.dat'            u ($1/1e14):5 w l lc 2 lt 1 lw 2 title 'Ymu',\
#'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):5 w l lc 2 lt 2 lw 2 title 'Ymu',\
#'APR_EOS_Cat.dat'            u ($1/1e14):6 w l lc 3 lt 1 lw 2 title 'Yn',\
#'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):6 w l lc 3 lt 2 lw 2 title 'Yn',\
#'APR_EOS_Cat.dat'            u ($1/1e14):7 w l lc 4 lt 1 lw 2 title 'Yp',\
#'SL_EOS_Cat+NV_smooth.dat'   u ($1/1e14):7 w l lc 4 lt 2 lw 2 title 'Yp',\



#2=dashed, 3=hashed, 4=dot, 5=dot-dash


