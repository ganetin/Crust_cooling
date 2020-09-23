#         nbar              rho           presion         ne            nu            nd            ns            mue           
#             1                 2             3             4             5             6             7 
plot 'strange_eos_03_140_150.dat' u 1:($6**(2./3.)-$5**(2./3.)-$4**(2./3.)) w l
replot 'strange_eos_03_140_150.dat' u 1:($5**(2./3.)-$7**(2./3.)-$4**(2./3.)) w l
