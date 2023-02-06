
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import mpl_scatter_density
from matplotlib.colors import LogNorm

map4=fits.open('integral_map_4.fits')
map16=fits.open('integral_map_16.fits')
map32=fits.open('integral_map_32.fits')
map64=fits.open('integral_map_64.fits')

inte4=map4[0].data
inte16=map16[0].data
inte32=map32[0].data
inte64=map64[0].data

#Select resolution
bool4=True
bool16=False
bool32=False
bool64=False


#4

if (bool4):
  plt.matshow(inte4,norm = LogNorm(vmin=np.min(inte4), vmax = 1),cmap = 'inferno',extent=[0,500,0,500])
  cbar=plt.colorbar()
  cbar.set_label(r'I [erg cm$^{-2} s^{-1} st^{-1}$]')

  plt.title('resolution= 4x4 pixels ')
  plt.xlabel('x-pixels [Mm]')
  plt.ylabel('y-pixels [Mm]')
  plt.savefig('4x4.eps',format='eps')


#16

if (bool16):
  plt.matshow(inte16,norm = LogNorm(vmin=np.min(inte16), vmax = 1),cmap = 'inferno',extent=[0,500,0,500])
  cbar=plt.colorbar()
  cbar.set_label(r'I [erg cm$^{-2} s^{-1} st^{-1}$]')

  plt.title('resolution= 16x16 pixels ')
  plt.xlabel('x-pixels [Mm]')
  plt.ylabel('y-pixels [Mm]')
  plt.savefig('16x16.eps',format='eps')


#32

if (bool32):
  plt.matshow(inte32,norm = LogNorm(vmin=np.min(inte32), vmax = 1),cmap = 'inferno',extent=[0,500,0,500])
  cbar=plt.colorbar()
  cbar.set_label(r'I [erg cm$^{-2} s^{-1} st^{-1}$]')

  plt.title('resolution= 32x32 pixels ')
  plt.xlabel('x-pixels [Mm]')
  plt.ylabel('y-pixels [Mm]')
  plt.savefig('32x32.eps',format='eps')


#64

if (bool64):
  plt.matshow(inte64,norm = LogNorm(vmin=np.min(inte64), vmax = 1),cmap = 'inferno',extent=[0,500,0,500])
  cbar=plt.colorbar()
  cbar.set_label(r'I [erg cm$^{-2} s^{-1} st^{-1}$]')

  plt.title('resolution= 64x64 pixels ')
  plt.xlabel('x-pixels [Mm]')
  plt.ylabel('y-pixels [Mm]')
  plt.savefig('64x64.eps',format='eps')



