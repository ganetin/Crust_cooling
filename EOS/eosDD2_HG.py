import numpy as np
from matplotlib import rc
from scipy import interpolate

Pr,rhor,nbr,nnr,npr,ner = np.loadtxt("DD2_HG.dat", usecols=(1,2,3,4,5,8), unpack=True)

density = np.linspace(0.01,max(nbr),200)

iP = interpolate.splrep(nbr,Pr)
P = interpolate.splev(density, iP)

irho = interpolate.splrep(nbr,rhor)
rho = interpolate.splev(density, irho)

inn = interpolate.splrep(nbr,nnr)
nn = interpolate.splev(density, inn)

inp = interpolate.splrep(nbr,npr)
np = interpolate.splev(density, inp)

ine = interpolate.splrep(nbr,ner)
ne = interpolate.splev(density, ine)

nb=density
P=P*1.6022e33
rho=rho*1.7827e12
Ye=ne/density
Yp=np/nb
Ym=Yp-Ye
Yn=1.-Yp

with open('eosDD2_HG.dat','w') as g:
  for i in range(len(nb)):
   if Ym[i]>0.:
      g.write(\
      "{0:16.6e}".format(float(rho[i]))+\
      "{0:16.6e}".format(float(P[i]))+\
      "{0:16.6e}".format(float(nb[i]))+\
      "{0:16.6e}".format(float(Ye[i]))+\
      "{0:16.6e}".format(float(Ym[i]))+\
      "{0:16.6e}".format(float(Yn[i]))+\
      "{0:16.6e}".format(float(Yp[i]))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(1.))+\
      "{0:16.6e}".format(float(1.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "\n")
   else:
      g.write(\
      "{0:16.6e}".format(float(rho[i]))+\
      "{0:16.6e}".format(float(P[i]))+\
      "{0:16.6e}".format(float(nb[i]))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(1.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(1.))+\
      "{0:16.6e}".format(float(1.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "{0:16.6e}".format(float(0.))+\
      "\n")
