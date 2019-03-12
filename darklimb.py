#-*- coding: utf-8
"""
	Program to correct limb darkening in python based on ILD/swwidl
	Calling sequence:
		import darkimb
		...
		...
	
	Coefficients taken from 
		Cox, A. N.: Allen's Astrophysical Quantities, Springer, 2000
		taken from IDL
	Author:	Angel Daniel Martínez-Cifuentes
	Date:	march, 10 2019
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os


#Define parameters
wavelenght = 6173.0 # Value in Angstroms
name = 'example.fits'

ll =1.0* wavelenght

data = fits.getdata(name,0)
head = fits.getheader(name,0)

xcen = head['CRPIX1']
ycen = head['CRPIX2']
radius = head['RSUN_OBS']/head['CDELT1'] #In pixels
size = head['NAXIS1']


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

xarr = np.arange(0,size,1.)
yarr = np.arange(0,size,1.)
xx, yy = np.meshgrid(xarr, yarr)#, sparse=True)
z = np.sqrt((xx-xcen)**2 + (yy-ycen)**2)
grid = z/radius
out = np.where(grid>1.0)
grid[out] = 0.0

#Equation
limbfilt =  1.0-ul-vl+ul*np.cos(np.arcsin(grid))+vl*np.cos(np.arcsin(grid))**2

imgout=np.array(array/limbfilt)

def writefits():
	os.system("rm -r limbcorrect.fits")
	hdu = fits.PrimaryHDU(imgout)
	hdul = fits.HDUList([hdu])
	hdul.writeto('limbcorrect.fits')

#writefits()

def figure():
	plt.figure(figsize=(8,8))
	plt.imshow(imgout,cmap = 'Greys_r',origin='lower')
	title = "corrected"
	plt.title(title)
	plt.savefig(title+".png")
	plt.show()



figure()

