# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 14:20:23 2017

@author: aouyed
"""
import matplotlib.pyplot as plt
from natsort import natsorted

from numpy import *
import numpy as np
import scipy as sp
from scipy import integrate
from scipy import interpolate

i=326
M=25
f=10**(1/M)
nref=1e-12
print(nref*f**(i-1))

j=81
Mt=25

nref=1.0000000e-01
f=10**(1/Mt)
print(nref*f**(j-1))