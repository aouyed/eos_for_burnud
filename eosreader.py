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

filename='eos2.ther,o'
eos=np.loadtxt(filename,skiprows=1)
nb=0.5
M=2
f=10**(1/M)
nref=1e-12
index=log(nb/nref)/log(f)+1
i=int(index)
print("i")
print(i)
nb1=nref*f**(i+1-1)
nb0=nref*f**(i-1)
y0=eos[i,3]
y1=eos[i+1,3]
y=y0+(nb-nb0)*(y1-y0)/(nb1-nb0)
print(y)
print(y1)
print(y0)
#mB=(eos[i,5]+1)*938
#print(eos[i,5])