
import matplotlib.pyplot as plt
import numpy as num
from scipy import interpolate
from scipy.interpolate import interp1d
import matplotlib.patches as patches
import os  


def read_Teff(filename,delay):   
    time=[]
    teff=[]
    g=open(filename,"r") #OK
    lines = g.readlines()[26:]
    for x in lines:
        time.append(float(x.split()[1])-12.5-delay)
        teff.append(float(x.split()[2])/1e6)
    g.close()    
    return time,teff

fig = plt.figure()

ax = fig.add_subplot(321)
ax.axis([10,15, 0,1.5])

#ax2.set_xlabel(r"$log10(\rho) [g cm$^{-3}$])", fontsize=13)
#ax.set_ylabel(r"Heat [MeV/nucleon]", fontsize=13)
plt.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=True,direction='in',left=True,right=True,labelright=False,labelleft=True)

ax.minorticks_on()

plt.tick_params(axis='x',which='both',bottom=True,top=True,labelbottom=False,direction='in')
plt.tick_params(axis='y',which='both',left=True,right=True,labelleft=True,direction='in')
for label in ax.get_yticklabels()+ax.get_xticklabels():
        label.set_fontsize(14)

rho,heat= num.loadtxt("heat_sources_HZ08.dat", usecols=(0,1), unpack=True)
for i in range(0,len(rho),2):
 ax.plot(num.log10(float(rho[i])),heat[i], 'ro')
 ax.vlines(num.log10(float(rho[i])), ymin=heat[i], ymax=heat[i+1],color='red')
ax.axvline(num.log10(1.4e14), color='black', linestyle=':', linewidth=1.5)
ax.annotate(r'HZ08',xy=(10.2,1.0),fontsize=16,color='black')

ax3 = fig.add_subplot(323)
ax3.axis([10,15, 0,1.5])

ax3.set_ylabel(r"Heat [MeV/nucleon]", fontsize=13)
plt.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=True,direction='in',left=True,right=True,labelright=False,labelleft=True)

ax3.minorticks_on()

plt.tick_params(axis='x',which='both',bottom=True,top=True,labelbottom=False,direction='in')
plt.tick_params(axis='y',which='both',left=True,right=True,labelleft=True,direction='in')
for label in ax3.get_yticklabels()+ax3.get_xticklabels():
        label.set_fontsize(14)

rho,heat= num.loadtxt("heat_sources_BSk20.dat", usecols=(0,1), unpack=True)
for i in range(0,len(rho),2):
 ax3.plot(num.log10(float(rho[i])),heat[i], 'ro')
 ax3.vlines(num.log10(float(rho[i])), ymin=heat[i], ymax=heat[i+1],color='red')
ax3.axvline(num.log10(1.4e14), color='black', linestyle=':', linewidth=1.5)
ax3.annotate(r'BSk20',xy=(10.2,1.0),fontsize=16,color='black')

ax3.annotate(r'crust',xy=(13.85,0.4),fontsize=10,color='black',rotation=90)
ax3.annotate(r'core',xy=(14.15,0.341),fontsize=10,color='black',rotation=90)

ax2 = fig.add_subplot(325)
# axis labels 
ax2.axis([10,15, 0,1.5])

ax2.set_xlabel(r"log10($\rho$ [g cm$^{-3}$])", fontsize=13)
#ax2.set_ylabel(r"Heat [MeV/nucleon]", fontsize=13)
plt.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=True,direction='in',left=True,right=True,labelright=False,labelleft=True)

ax2.minorticks_on()


plt.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=True,direction='in',left=True,right=True,labelright=False,labelleft=True)
rho,heat= num.loadtxt("heat_sources_BSk21.dat", usecols=(0,1), unpack=True)
for i in range(0,len(rho),2):
 ax2.plot(num.log10(float(rho[i])),heat[i], 'ro')
 ax2.vlines(num.log10(float(rho[i])), ymin=heat[i], ymax=heat[i+1],color='red')
ax2.axvline(num.log10(1.4e14), color='black', linestyle=':', linewidth=1.5)
ax2.annotate(r'BSk21',xy=(10.2,1.0),fontsize=16,color='black')


for label in ax2.get_yticklabels()+ax2.get_xticklabels():
        label.set_fontsize(14)


ax2.minorticks_on()


ax4 = fig.add_subplot(122)
# axis labels 
ax4.axis([-1,20, 0.6,1.3])

ax4.set_xlabel(r"time since accretion ended [year]", fontsize=11)
ax4.set_ylabel(r"$T^\infty$ [MK]", fontsize=13)
plt.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=True,direction='in',left=True,right=True,labelright=True,labelleft=False)

ax4.minorticks_on()
ax4.yaxis.set_label_position("right")

for label in ax4.get_yticklabels()+ax4.get_xticklabels():
        label.set_fontsize(14)
ax4.minorticks_on()


plt.tick_params(axis='both',which='both',bottom=True,top=True,labelbottom=True,direction='in',left=True,right=True,labelright=True,labelleft=False)


time,teff,errteff= num.loadtxt("../Observations/KS_1731-260.dat", usecols=(0,1,2), unpack=True)
ax4.errorbar((time-51930.5)/365.25,(teff/8.617343e-5)/1e6, yerr=errteff*2/8.617343e-5/1e6, fmt='ko', capsize=5)


colors=['grey','brown','orange','olive','green','cyan','blue','purple','pink','red']
styles=['-','--',':','-.']

Tempok=[]
Mdotok=[]
Chiok=[]
delay=5e3
j=-1
temp=5.e7
temp=temp-.1e7
while temp<=5.1e7:
 temp=temp+.1e7
 Temp="{:.1e}".format(temp) 
 j=j+1
 i=-1
 mdot=2.0e-9
 mdot=mdot-0.1e-9
 while mdot<=2.1e-9:
  mdot=mdot+0.1e-9
  Mdot="{:.1e}".format(mdot) 
  i=i+1
  name="Teff_"+Temp+"_"+Mdot+"_0_Acc_5e3_1e8_1.4"
  print(name)
  if os.path.isfile(name+".dat")==True:
   print(Temp,Mdot)

   f = interpolate.interp1d(read_Teff(name+".dat",delay)[0],read_Teff(name+".dat",delay)[1]) 
   T0=f((time[0]-51930.5)/365.25)
   T0max=(teff[0]+errteff[0]*2)/8.617343e-5/1e6
   T0min=(teff[0]-errteff[0]*2)/8.617343e-5/1e6
   T7=f((time[7]-51930.5)/365.25)
   T7max=(teff[7]+errteff[7]*2)/8.617343e-5/1e6
   T7min=(teff[7]-errteff[7]*2)/8.617343e-5/1e6
   chi2=0
   for k in range(len(time)):
      Ti=(f((time[k]-51930.5)/365.25))
      Tiobs=(teff[k])/8.617343e-5/1e6
      Timax=(teff[k]+errteff[k]*2)/8.617343e-5/1e6
      Timin=(teff[k]-errteff[k]*2)/8.617343e-5/1e6
      erri=errteff[k]*2/8.617343e-5/1e6
      chi2=chi2+((Ti-Tiobs)/erri)**2

   Tempok.append(temp)
   Mdotok.append(mdot)
   Chiok.append(chi2**0.5)
   print("**",Temp,Mdot,chi2**0.5)

index_optimal=Chiok.index(min(Chiok))
print("Optimal model:",Tempok[index_optimal],Mdotok[index_optimal],Chiok[index_optimal])
name="Teff_"+"{:.1e}".format(Tempok[index_optimal])+"_"+"{:.1e}".format(Mdotok[index_optimal])+"_0_Acc_5e3_1e8_1.4"
ax4.plot(read_Teff(name+".dat",delay)[0],read_Teff(name+".dat",delay)[1],'-',color='k',linestyle='-',linewidth=1.5,label="{:.1e}".format(Tempok[index_optimal])+"_"+"{:.1e}".format(Mdotok[index_optimal]))



l1=ax4.legend(loc='upper right',fontsize=11,ncol=1,borderaxespad=0.1, labelspacing=0.1,
handletextpad=+0.18, markerscale=.5,columnspacing=0.5,handlelength=1.7)
l1.draw_frame(False)


plt.savefig("DCH_KS.pdf", format="pdf", bbox_inches="tight")


