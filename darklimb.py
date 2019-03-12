#-*- coding: utf-8
"""
	Program to correct limb darkening in python based on ILD/swwidl
	Calling sequence:
		import darkimb
		...
		...
	Required libraries: 
		Numpy
		Astropy
		Matplotlib
		If not: pip install <library>
	Author:	Angel Daniel Martínez-Cifuentes
	Date:	march, 10 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math as m
#Define parameters

wavelenght = 6173.0 # Value in Angstroms
name = 'example.fits'


ll =1.0* wavelenght

data = fits.getdata(name,0)
head = fits.getheader(name,0)

#Define darkilmb_u

def darklimb_u(ll):
	pll = np.array([1.0,ll,ll**2,ll**3,ll**4,ll**5])
	au = -8.9829751
	bu = 0.0069093916
	cu = -1.8144591e-6
	du = 2.2540875e-10
	eu = -1.3389747e-14
	fu = 3.0453572e-19
	a=np.array([au,bu,cu,du,eu,fu])
	ul = sum(a*pll)
	return ul

def darklimb_v(ll):
	pll = np.array([1.0,ll,ll**2,ll**3,ll**4,ll**5])
	av = 9.2891180
	bv = -0.0062212632
	cv = 1.5788029e-6
	dv = -1.9359644e-10
	ev = 1.1444469e-14
	fv = -2.599494e-19
	a=np.array([av,bv,cv,dv,ev,fv])
	vl = sum(a*pll)
	return vl

array = np.array(data)
NaNs = np.isnan(array)
array[NaNs] = 0.0

ul = darklimb_u(ll)
vl = darklimb_v(ll)

dist = np.loadtxt('d.dat')
print(type(dist))




exit()
limbfilt =  1.0-ul-vl+ul*np.cos(np.arcsin(dist_grid))+vl*np.cos(np.arcsin(dist_grid))**2
imgout=np.array(newdata/limbfilt)


plt.imshow(imgout,cmap = 'Greys_r',origin='lower')
plt.show()





















