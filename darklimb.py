#-*- coding: utf-8
"""
	Program to correct limb darkening in python based on ILD/swwidl
	Based on:
	 https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/darklimb_correct.pro

	Coefficients taken from 
		Cox, A. N.: Allen's Astrophysical Quantities, Springer, 2000
		taken from IDL
	Author:	Angel 
	Date:	march, 10 2019
"""


import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os

def writefits(image):
	os.system("rm -r limbcorrect.fits")
	hdu = fits.PrimaryHDU(image)
	hdul = fits.HDUList([hdu])
	hdul.writeto('limbcorrect.fits')

def figure(image,title="image",save=False):
	plt.figure(figsize=(8,8))
	plt.title('Image')
	plt.imshow(image,cmap = 'Greys_r',origin='lower')
	if save == True:
		plt.savefig(title+".png")
	else:
		pass
	plt.show()


def darklimb(array):
	"""
	  Darklimb function
	  
	  Output:
	    Dos arrays: El primero el corregido, el segundo el original
	"""
	
	# Se hacen dos funciones, darklimb_u y darklimb_v, para aplicar
	# posteriormente el corregimiento 
	
	# -------------------------------------------------------------
	
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
	
	# -------------------------------------------------------------
	
	# Leer los archivos .fits (tabla y header)
	data = fits.getdata(name,0)
	head = fits.getheader(name,0)

	# Parametros necesarios para la función
	wavelength = head['WAVELNTH'] # Longitud de onda
	xcen = head['CRPIX1'] # Centro de la imagen en x
	ycen = head['CRPIX2'] # Centro de la imagen en y
	radius = head['RSUN_OBS']/head['CDELT1'] # Resultado en Pixels
	size = head['NAXIS1'] # Tamaño del array en x

	ll =1.0*wavelength

	array = np.array(data) # Convertir los datos a arrays de numpy
	NaNs = np.isnan(array) # Mirar donde los datos son NAN
	array[NaNs] = 0.0      # Hacer cero todos los NAN

	# Aplicar los corregimientos
	ul = darklimb_u(ll)
	vl = darklimb_v(ll)

	xarr = np.arange(0,size,1.0) # Hacer array en x
	yarr = np.arange(0,size,1.0) # Hacer array en y
	xx, yy = np.meshgrid(xarr, yarr) # Hacer array en xy (cuadrado)
	# z: Hacer array cuyo centro en cero es el centro xy del array
	# Se hace a modo de realizar un circulo 
	z = np.sqrt((xx-xcen)**2 + (yy-ycen)**2) 
	# grid: Normalizar el circulo tal que adentro sea menor a 1
	grid = z/radius 
	out = np.where(grid>1.0) # Buscar donde grid es mayor a 1
	grid[out] = 0.0 # Hacer cero donde sea mayor a 1

	# Equation: Hacer ecuación para aplicar el corregimiento del limbo
	# La ecuación se toma a modo de corregir limbo
	limbfilt =  1.0-ul-vl+ul*np.cos(np.arcsin(grid))+vl*np.cos(np.arcsin(grid))**2

	# Imagen final
	imgout=np.array(array/limbfilt)
	
	return imgout, array


name = 'hmi.ic_45s.2014.02.04_03_44_15_TAI.continuum.fits'
corrected, original = darklimb(name)

#writefits(corrected)
figure(corrected,title='corrected1',save=True)
figure(original,title='original1',save=True)

