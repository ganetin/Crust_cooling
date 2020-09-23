import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from matplotlib.ticker import MultipleLocator
import matplotlib
matplotlib.colors.ColorConverter.colors['z'] = (0.4,0.4,0.4)
from scipy import interpolate
from matplotlib import colors
from scipy.optimize import curve_fit
def func(x, a, b, c):
    return a*np.exp(-x/b) + c
def func2(x, a, b):
    return a*(x)**b


t1,Teff1,dTeff1 = np.loadtxt('KS_1731-260.dat'    , usecols=(0,1,2), unpack=True)
t01=51930.5
t2,Teff2,dTeff2 = np.loadtxt('MXB_1659-29.dat'    , usecols=(0,1,2), unpack=True)
t02=52159.5
t3,Teff3,dTeff3= np.loadtxt('EXO_0748-676_2.dat'    , usecols=(1,2,3), unpack=True)
t03=54714 
t4,Teff4,dTeff4= np.loadtxt('XTE_J1701-462_2.dat'    , usecols=(1,2,3), unpack=True)
t4b,Teff4b,dTeff4b= np.loadtxt('XTE_J1701-462_2b.dat'    , usecols=(1,2,3), unpack=True)
t04= 54321.95
t5,Teff5,dTeff5= np.loadtxt('IGR_J17480.dat'    , usecols=(0,1,2), unpack=True)
t05=55556
t6,Teff6,dTeff6= np.loadtxt('Swift_J17480.dat'    , usecols=(0,1,2), unpack=True)
t06=56166
t7,Teff7,dTeff7= np.loadtxt('MAXI_J0556-332.dat'    , usecols=(1,2,3), unpack=True)
t7b,Teff7b,dTeff7b= np.loadtxt('MAXI_J0556-332b.dat'    , usecols=(1,2,3), unpack=True)


## PLOTS
rc('font', family='serif')
rc(('ytick.major','ytick.minor'), pad=7)
fig = plt.figure()
ax = fig.add_subplot(121)
for label in ax.get_yticklabels():
        label.set_fontsize(18)
for label in ax.get_xticklabels():
        label.set_fontsize(18)
# Second plot: L
ax.axis([1,6000,0.2,2])
ax.set_autoscale_on(False)
ax.xaxis.labelpad = 5
ax.yaxis.labelpad = 5
ax.axes.get_xaxis().set_visible(True)
ax.get_yticks(minor=True)
ax.minorticks_on()
fig.text(0.5, 0.,r"Time since last outburst (days)" ,fontsize=22, ha='center')
#ax.set_xlabel("")
ax.set_ylabel(r"$T_{\rm eff}^\infty$ (MK)", fontsize=22)
ax.set_xscale('log')
x=np.logspace(0,4)

for i in range(len(t1)):
 ax.errorbar((t1[i]-t01),((Teff1[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff1[i]/8.617343e-5/1e6]) ,markeredgecolor='k',markerfacecolor = 'k',mew=1.5,markersize=1,ecolor='k')
for i in range(1):
 ax.errorbar((t1[i]-t01),((Teff1[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff1[i]/8.617343e-5/1e6]) ,markeredgecolor='k',markerfacecolor = 'k',mew=1.5,markersize=1,ecolor='k',label="KS")# - Merritt+16")
popt, pcov = curve_fit(func, (t1-t01),Teff1,sigma=dTeff1,p0=[500,100,40])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'k')
popt2, pcov2 = curve_fit(func2, (t1-t01),Teff1,sigma=dTeff1)
print("KS &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))


for i in range(len(t2)):
 ax.errorbar((t2[i]-t02),((Teff2[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff2[i]/8.617343e-5/1e6]) ,markeredgecolor='r',markerfacecolor = 'r',mew=1.5,markersize=1,ecolor='r')
for i in range(1):
 ax.errorbar((t2[i]-t02),((Teff2[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff2[i]/8.617343e-5/1e6]) ,markeredgecolor='r',markerfacecolor = 'r',mew=1.5,markersize=1,ecolor='r',label="MXB")# - Cackett+13")
popt, pcov = curve_fit(func, (t2-t02),Teff2,sigma=dTeff2,p0=[500,100,40])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'r')
popt2, pcov2 = curve_fit(func2, (t2-t02),Teff2,sigma=dTeff2)
print("MXB &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))



for i in range(len(t3)):
 ax.errorbar((t3[i]-t03),((Teff3[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff3[i]/8.617343e-5/1e6]) ,markeredgecolor='b',markerfacecolor = 'b',mew=1.5,markersize=1,ecolor='b')
for i in range(1):
 ax.errorbar((t3[i]-t03),((Teff3[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff3[i]/8.617343e-5/1e6]) ,markeredgecolor='b',markerfacecolor = 'b',mew=1.5,markersize=1,ecolor='b',label="EXO")# - Degenaar+14")
ax.axhline(y=(94.6+5.6)/8.617343e-5/1e6,ls='--', linewidth=1, color = 'b')
ax.axhline(y=(94.6-16.0)/8.617343e-5/1e6,ls='--', linewidth=1, color = 'b')
popt, pcov = curve_fit(func, (t3-t03),Teff3,sigma=dTeff3,p0=[500,100,40])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'b')
popt2, pcov2 = curve_fit(func2, (t3-t03),Teff3,sigma=dTeff3)
print("EXO &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))


for i in range(len(t4)):
 ax.errorbar((t4[i]),((Teff4[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff4[i]/8.617343e-5/1e6]) ,markeredgecolor='orange',markerfacecolor = 'orange',mew=1.5,markersize=1,ecolor='orange')
for i in range(1):
 ax.errorbar((t4[i]),((Teff4[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff4[i]/8.617343e-5/1e6]) ,markeredgecolor='orange',markerfacecolor = 'orange',mew=1.5,markersize=1,ecolor='orange',label="XTE")# - Fridriksson+11")
popt, pcov = curve_fit(func,  (t4b),Teff4b,sigma=dTeff4b,p0=[500,100,40])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'orange')
popt2, pcov2 = curve_fit(func2, (t4b),Teff4b,sigma=dTeff4b)
print("XTE &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))

for i in range(len(t5)):
 ax.errorbar((t5[i]-t05),((Teff5[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff5[i]/8.617343e-5/1e6]) ,markeredgecolor='green',markerfacecolor = 'green',mew=1.5,markersize=1,ecolor='green')
for i in range(1):
 ax.errorbar((t5[i]-t05),((Teff5[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff5[i]/8.617343e-5/1e6]) ,markeredgecolor='green',markerfacecolor = 'green',mew=1.5,markersize=1,ecolor='green',label="IGR")# J17480-2446 - Degenaar+15")
ax.axhline(y=(73.6+1.6)/8.617343e-5/1e6,ls='--', linewidth=1, color = 'green')
ax.axhline(y=(73.6-1.6)/8.617343e-5/1e6,ls='--', linewidth=1, color = 'green')
popt, pcov = curve_fit(func, (t5-t05),Teff5,sigma=dTeff5,p0=[200,30,120])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'green')
popt2, pcov2 = curve_fit(func2,(t5-t05),Teff5,sigma=dTeff5)
print("IGR &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))

l1=ax.legend(loc='lower left',numpoints=1,fontsize=12,ncol=2)#,borderaxespad=0.1, numpoints=1,labelspacing=0.1,handletextpad=+0.18, markerscale=1.,ncol=1,fontsize=15)
l1.draw_frame(False)


#plt.savefig("Thermal_relax.eps", format="eps", bbox_inches="tight")




### PLOTS
#rc('font', family='serif')
#rc(('ytick.major','ytick.minor'), pad=7)
#fig = plt.figure()
ax = fig.add_subplot(122)
for label in ax.get_yticklabels():
        label.set_fontsize(18)
for label in ax.get_xticklabels():
        label.set_fontsize(18)
# Second plot: L
ax.axis([1,1000,0.2,4.5])
ax.set_autoscale_on(False)
ax.xaxis.labelpad = 5
ax.yaxis.labelpad = 5
ax.axes.get_xaxis().set_visible(True)
ax.get_yticks(minor=True)
ax.minorticks_on()
#ax.set_xlabel(r"Time since last outburst (days)" ,fontsize=22)
#ax.set_ylabel(r"$T_{\rm eff}^\infty$ (MK)", fontsize=22)
ax.set_xscale('log')


for i in range(len(t6)):
 ax.errorbar((t6[i]-t06),((Teff6[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff6[i]/8.617343e-5/1e6]) ,markeredgecolor='k',markerfacecolor = 'k',mew=1.5,markersize=1,ecolor='k')
for i in range(1):
 ax.errorbar((t6[i]-t06),((Teff6[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff6[i]/8.617343e-5/1e6]) ,markeredgecolor='k',markerfacecolor = 'k',mew=1.5,markersize=1,ecolor='k',label="Swift")# J174805 - Degenaar+15")
ax.axhline(y=(89.7+1.7)/8.617343e-5/1e6,ls='--', linewidth=1, color = 'k')
ax.axhline(y=(89.7-1.7)/8.617343e-5/1e6,ls='--', linewidth=1, color = 'k')
popt, pcov = curve_fit(func, (t6-t06),Teff6,sigma=dTeff6,p0=[200,30,120])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'k')
popt2, pcov2 = curve_fit(func2,(t6-t06),Teff6,sigma=dTeff6)
print("Swift &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))


for i in range(len(t7)):
 ax.errorbar((t7[i]),((Teff7[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff7[i]/8.617343e-5/1e6]) ,markeredgecolor='b',markerfacecolor = 'b',mew=1.5,markersize=1,ecolor='b')
for i in range(1):
 ax.errorbar((t7[i]),((Teff7[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff7[i]/8.617343e-5/1e6]) ,markeredgecolor='b',markerfacecolor = 'b',mew=1.5,markersize=1,ecolor='b',label="MAXI")# J0556-332 - Homan+14")
popt, pcov = curve_fit(func, t7b,Teff7b,sigma=dTeff7b,p0=[200,30,120])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'b')
popt2, pcov2 = curve_fit(func2, t7b,Teff7b,sigma=dTeff7b)
print("MAXI &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))
#print Teff7,dTeff7b
t8,Teff8,dTeff8= np.loadtxt('Aql_X-1_a.dat'    , usecols=(0,1,2), unpack=True)
t08= 55925
for i in range(len(t8)):
 ax.errorbar((t8[i]-t08),((Teff8[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff8[i]/8.617343e-5/1e6]) ,markeredgecolor='purple',markerfacecolor = 'purple',mew=1.5,markersize=1,ecolor='purple')
for i in range(1):
 ax.errorbar((t8[i]-t08),((Teff8[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff8[i]/8.617343e-5/1e6]) ,markeredgecolor='purple',markerfacecolor = 'purple',mew=1.5,markersize=1,ecolor='purple',label="Aql")# X-1 - Waterhouse+16")
t8,Teff8,dTeff8= np.loadtxt('Aql_X-1_a2.dat'    , usecols=(0,1,2), unpack=True)
t08= 55925
popt, pcov = curve_fit(func, (t8-t08),Teff8,sigma=dTeff8,p0=[200,30,120])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'purple')
popt2, pcov2 = curve_fit(func2,(t8-t08),Teff8,sigma=dTeff8)
print("Aql X-1 a &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))
#popt, pcov = curve_fit(func, (t8-t08),Teff8,sigma=dTeff8,p0=[200,30,120])
#print 
#print
#print 'Aql X-1 - a'
#print '------------------'
#print 'Exponential decay'
#print 'tau=',"{:.2f}".format(popt[1]),'+/-',"{:.2f}".format(np.sqrt(np.diag(pcov))[1])
#print 'A=',"{:.2f}".format(popt[0]),'+/-',"{:.2f}".format(np.sqrt(np.diag(pcov))[0])
#print 'B=',"{:.2f}".format(popt[2]),'+/-',"{:.2f}".format(np.sqrt(np.diag(pcov))[2])
#ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'purple')
#print '------------------'
#print 'Power law'
#popt, pcov = curve_fit(func2, (t8-t08),Teff8,sigma=dTeff8)
#print 'alpha=',"{:.2f}".format(popt[0]),'+/-',"{:.2f}".format(np.sqrt(np.diag(pcov))[0])
#print 'beta=',"{:.3f}".format(popt[1]),'+/-',"{:.3f}".format(np.sqrt(np.diag(pcov))[1])

t8,Teff8,dTeff8= np.loadtxt('Aql_X-1_b.dat'    , usecols=(0,1,2), unpack=True)
t08= 56518 
for i in range(len(t8)):
 ax.errorbar((t8[i]-t08),((Teff8[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff8[i]/8.617343e-5/1e6]) ,markeredgecolor='plum',markerfacecolor = 'plum',mew=1.5,markersize=1,ecolor='plum')
popt, pcov = curve_fit(func, (t8-t08),Teff8,sigma=dTeff8,p0=[200,30,120])
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'plum')
popt2, pcov2 = curve_fit(func2,(t8-t08),Teff8,sigma=dTeff8)
print("Aql X-1 b &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))


t8,Teff8,dTeff8= np.loadtxt('Aql_X-1_c.dat'    , usecols=(0,1,2), unpack=True)
t08= 57100 
for i in range(len(t8)):
 ax.errorbar((t8[i]-t08),((Teff8[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff8[i]/8.617343e-5/1e6]) ,markeredgecolor='orchid',markerfacecolor = 'orchid',mew=1.5,markersize=1,ecolor='orchid')
#for i in range(1):
# ax.errorbar((t8[i]-t08),((Teff8[i]/8.617343e-5)/1e6), fmt='o',yerr=np.array([-dTeff8[i]/8.617343e-5/1e6]) ,markeredgecolor='orchid',markerfacecolor = 'orchid',mew=1.5,markersize=1,ecolor='orchid',label="Aql X-1 - Waterhouse+16")


t9,t9m,t9p,Teff9,Teff9m,Teff9p= np.loadtxt('1RXS.dat'    , usecols=(0,1,2,3,4,5), unpack=True)
ax.errorbar(t9, Teff9/8.617343e-5/1e6, xerr=[-t9m,t9p],yerr=[-Teff9m/8.617343e-5/1e6,Teff9p/8.617343e-5/1e6] , fmt='o',markeredgecolor='green',markerfacecolor = 'green',mew=1.5,markersize=1,ecolor='green')
for i in range(1):
 ax.errorbar(t9[i], Teff9[i]/8.617343e-5/1e6,yerr=np.array(-Teff9m[i]/8.617343e-5/1e6,Teff9p[i]/8.617343e-5/1e6), fmt='o',markeredgecolor='green',markerfacecolor = 'green',mew=1.5,markersize=1,ecolor='green',label="1RXS")# J1804 - Parikh+17")
popt, pcov = curve_fit(func,t9,Teff9,sigma=(-Teff9m+Teff9p)/2.,p0=[200,30,120])
#print Teff9,(-Teff9m+Teff9p)/2.
ax.plot(x,func(x,popt[0],popt[1],popt[2])/8.617343e-5/1e6, color = 'green')
popt2, pcov2 = curve_fit(func2, t9,Teff9,sigma=(-Teff9m+Teff9p)/2.)
print("1RXS &","{:.2f}".format(popt[1]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[1]),"& {:.2f}".format(popt[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[0]),"& {:.2f}".format(popt[2]),"& {:.2f}".format(np.sqrt(np.diag(pcov))[2]),"& & {:.2f}".format(popt2[0]),"& {:.2f}".format(np.sqrt(np.diag(pcov2))[0]),"& {:.3f}".format(popt2[1]),"& {:.3f}\\".format(np.sqrt(np.diag(pcov2))[1]))


l1=ax.legend(loc='lower left',numpoints=1,fontsize=12,ncol=2)#,borderaxespad=0.1, numpoints=1,labelspacing=0.1,handletextpad=+0.18, markerscale=1.,ncol=1,fontsize=15)
l1.draw_frame(False)



plt.savefig("Thermal_relax.eps", format="eps", bbox_inches="tight", orientation = 'portrait')




