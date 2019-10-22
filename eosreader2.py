# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 15:42:54 2017

@author: aouyed
"""

import matplotlib.pyplot as plt
from natsort import natsorted

from numpy import *
import numpy as np
import scipy as sp
from scipy import integrate
from scipy import interpolate

filename='eos2H.thermo'
mn=939.565413
fm3=(197.3269788)**3
fz=16
bot=0.15
rg=0.85

def quarkEOS(nB,T1):
    fm3=(197.3269788)**3
    mn=939.565413
    
    c = 3*10**10
    nBAsh = fm3*nB
    ne2=0.2*nBAsh;
    a_s=0.4
    mufactor=1.0 / (1.0 - 2.0 * a_s /np.pi)
    nU2 = ne2 + nBAsh
    nS2 = (0.5*nBAsh)
    MuS2 = (nS2*np.pi**2*mufactor)**(1/3)
    MuS2=np.sqrt(MuS2*MuS2+150**2)
    nD2 = 3*nBAsh - nS2 - nU2
    MuD2 = (mufactor*nD2*np.pi**2)**(1/3)
    MuU2 = (mufactor*nU2*np.pi**2)**(1/3)
    MuE2=0.0
    B=145**4
    p = (MuS2**4 + MuD2**4 + MuU2**4 + MuE2**4/3.0) + 2*np.pi**2 *T1**2 *(MuS2**2 + MuD2**2 + MuU2**2 + MuE2**2/3.0) + 205.641*T1**4
    p=p/(4*np.pi**2)-B
    p=p- ((7.0 / 60.0 * T1 * T1 * T1 * T1 * np.pi*np.pi * 50 * a_s / (21.0 * np.pi)) + 2.0 * a_s / np.pi * (0.5 * T1 * T1 * MuS2 * MuS2 + MuS2**4/ (4.0 * np.pi**2)));    
    p=p- ((7.0 / 60.0 * T1 * T1 * T1 * T1 * np.pi*np.pi * 50 * a_s / (21.0 * np.pi)) + 2.0 * a_s / np.pi * (0.5 * T1 * T1 * MuD2 * MuD2 + MuD2**4/ (4.0 * np.pi**2)));    
    p=p- ((7.0 / 60.0 * T1 * T1 * T1 * T1 * np.pi*np.pi * 50 * a_s / (21.0 * np.pi)) + 2.0 * a_s / np.pi * (0.5 * T1 * T1 * MuU2 * MuU2 + MuU2**4/ (4.0 * np.pi**2)));    

    u=3*p + 4*B;
    s=4*np.pi**2 *T1*(MuS2**2 + MuD2**2 + MuU2**2 + MuE2**2/3.0) + 4*205.641*T1**3;
    s=s/(4*np.pi**2) 
    
    f=u - s*T1
    #return(f/(nBAsh*mn)-1)
    return(f/(nB*fm3))
def eosReader(nBb, T, q):
    Nn=326
    M=25
    
    f=10**(1/M)
    nref=1e-12
    Tref=1.0000000e-01
    f=10**(1/M)
    jndex=log(T/Tref)/log(f)+1
    j=int((jndex))

    y1=Tref*f**(j+1-1)
    y0=Tref*f**(j-1)
    y2=T
    
    
    M=25
    f=10**(1/M)
    nref=1e-12
    index=log(nBb/nref)/log(f)+1
    i=int((index))

    x1=nref*f**(i+1-1)
    x0=nref*f**(i-1)
    x2=nBb
    
    
    z0=eos[(j)*Nn+i,q]
    z1=eos[(j)*Nn+i+1,q]
    z2=eos[(j-1)*Nn+i,q]
    z3=eos[(j-1)*Nn+i+1,q]
    
    Na=(x1-x2)*(y2-y0)/((x1-x0)*(y1-y0))
    Nb=(x2-x0)*(y2-y0)/((x1-x0)*(y1-y0))
    Nc=(x1-x2)*(y1-y2)/((x1-x0)*(y1-y0))
    Nd=(x2-x0)*(y1-y2)/((x1-x0)*(y1-y0))
    
    z4=z0*Na+z1*Nb+z2*Nc+z3*Nd
    
    return((z4+1)*mn)
    
    

eos=np.loadtxt(filename,skiprows=1)
print(eos.shape)

q=6
T=1
nB=0.4
print(quarkEOS(nB,T))
print(quarkEOS(nB,T)/eosReader(nB,T,q+2))
nbArray=np.arange(0.01,1,0.01)
TArray=np.arange(1,60,1)

print(TArray)

eosFunction= lambda x: eosReader(x, T, q+2)
eosFunctionq= lambda x: quarkEOS(x, T)

eosFunctionT= lambda x: eosReader(nB, x, q+2)
eosFunctionqT= lambda x: quarkEOS(nB, x)
vecfunc=np.vectorize(eosFunction)
vecfuncq=np.vectorize(eosFunctionq)
vecfuncT=np.vectorize(eosFunctionT)
vecfuncqT=np.vectorize(eosFunctionqT)
qArray=vecfunc(nbArray)
qArrayq=vecfuncq(nbArray)
qArrayT=vecfuncT(TArray)
qArrayqT=vecfuncqT(TArray)


plt.figure()

plt.subplots_adjust(bottom=bot,right=rg)
plt.plot(TArray,qArrayT,linestyle='-',label=r'$\rm{hadronic}$',color='blue')
plt.plot(TArray,qArrayqT,linestyle='--',label=r'$\rm{(u,d,s)}$',color='black')
plt.xlabel(r'$\rm{baryon}\;\rm{Temperature}\;[\rm{MeV}]$',fontsize=fz)
plt.ylabel(r'$\rm{free}\;\rm{energy}\;\rm{per}\;\rm{baryon}\;[\rm{MeV}]$',fontsize=fz)
plt.xticks(fontsize=fz);
plt.yticks(fontsize=fz);
plt.legend(loc='lower left',fontsize=fz)
tabA=np.array([TArray,np.array(qArrayT)])



#print(nbArray)


#np.savetxt('tableAmir.txt',np.c_[TArray,qArrayT])

#mB=(eos[i,5]+1)*938
#print(eos[i,5])