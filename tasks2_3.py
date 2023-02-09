import os
os.environ['XUVTOP']='C:\\Users\\sarit\\CHIANTI'
os.environ["HOME"]="C:\\Users\\sarit\\CHIANTI"


import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.stats import kde
import mpl_scatter_density

#plt.close('all')

#Functions
def find_closest(arr, val):  # Finding the nearest point in an array
    diff = abs(arr - val)
    closest = np.min(diff)
    index = np.where(diff == closest)[0][0]
    return index

#Importing data from .fits
hdulist1=fits.open('lgtg_400.fits')
hdulist2=fits.open('lgne_400.fits')
temp=hdulist1[0].data
z=hdulist1[1].data
nel=hdulist2[0].data #MKS units (ej. m−3)


#displaying the .fits
plt.figure(1)
image= plt.imshow(temp[:, :, 302],origin='lower', cmap = 'inferno')
cbar = plt.colorbar(image)
cbar.set_label('logT [cuentas]')
#plt.title(' ')
plt.xlabel('x-pixels [Mm]')
plt.ylabel('y-pixels [Mm]')


#------Task 2------#
#First acquaintance with the results of a numerical 3D experiment which are separately provided via data boxes.


#Booleans to execute for temperature or e- density and task 3
booltemp=False
boolne=False
booltask3=False


#Temperature

if (booltemp):
  plt.figure(2)

  #Mean temperature for z
  Tmean_h=np.zeros(len(z))
  for i in range(len(z)):
      Tmean_h[i]=np.mean(temp[i,:,:])

  Tmean_h=10**(Tmean_h)
  plt.plot(z,Tmean_h,label='T horizontal media')
  plt.yscale("log") 
  plt.xlabel('Heigth [Mm]')
  plt.ylabel('T [K] ')
  plt.xlim(-2.9,15)


  t_sec=10**temp[:,:,:]
  z_sec=z
  t_secf=t_sec.flatten()
  zf=[]
  for k in range(len(z_sec)):
      for i in range(len(t_sec[0,:,0])):
          for j in range(len(t_sec[0,0,:])):
              zf.append(z_sec[k])

  #Normalizing the data
  from astropy.visualization import LogStretch
  from astropy.visualization.mpl_normalize import ImageNormalize
  norm = ImageNormalize(vmin=0., vmax=10000, stretch=LogStretch())

  #Scatter plot of T vs height  
  fig3 = plt.figure(3)
  ax = fig3.add_subplot(1, 1, 1, projection='scatter_density')
  ax.scatter_density(zf,t_secf ,color='black',norm=norm)
  plt.plot(z,Tmean_h,"r-", lw=1,label='T horizontal media')
  plt.yscale("log") 
  plt.xlabel('Heigth [Mm]')
  plt.ylabel('T [K] ')
  plt.xlim(-2.9,15)
  plt.ylim(1e3,1e7)
  plt.legend()

  #plt.savefig('task2a_2.png',format='png')

  
#Electron density
if (boolne):
  plt.figure(4)
  
  #Mean electron density for z
  Nmean_h=np.zeros(len(z))
  for i in range(len(z)):
      Nmean_h[i]=np.mean(nel[i,:,:])


  Nmean_h=10**(Nmean_h)
  plt.plot(z,Nmean_h,label='N horizontal media')
  plt.yscale("log") 
  plt.xlabel('Heigth [Mm]')
  plt.ylabel(r'Ne [$ m^{−3}$]')
  plt.xlim(-2.9,15)
  #plt.ylim(1000,10000000)



  n_sec=10**nel[:,:,:]
  z_sec=z
  n_secf=n_sec.flatten()
  zf=[]
  for k in range(len(z_sec)):
      for i in range(len(n_sec[0,:,0])):
          for j in range(len(n_sec[0,0,:])):
              zf.append(z_sec[k])

  #Normalizing the data
  from astropy.visualization import LogStretch
  from astropy.visualization.mpl_normalize import ImageNormalize
  norm = ImageNormalize(vmin=0., vmax=10000, stretch=LogStretch())
  
  #Scatter plot of Ne vs height
  fig3 = plt.figure(5)
  ax = fig3.add_subplot(1, 1, 1, projection='scatter_density')
  ax.scatter_density(zf,n_secf ,color='black',norm=norm)
  plt.plot(z,Nmean_h,"r-", lw=1,label='Ne horizontal media')
  plt.yscale("log") 
  plt.xlabel('Heigth [Mm]')
  plt.ylabel(r'Ne [$ m^{−3}$]')
  plt.xlim(-2.9,15)
  #plt.ylim(1000,10000000)
  plt.legend()
  #plt.savefig('task2b_1.png',format='png')


#--------Task 3--------#
#we carry out the synthesis of the actual intensities emitted along the line of sigh using Chianti


if (booltask3):
  
  #To simplify the calculation: heights greater than 2 Mm above the photosphere
  z_sec=z[z>2]

  #Defining the resolution of the results
  resolution=4 
  a=np.linspace(0, 503, resolution )
  a=a.astype(int)


  T_sec3=10.**temp[180:,:,:]
  n_sec3=(10.**nel[180:,:,:])*1.e-6 #cm-3

  integral_map=[]
  for i in a:
      for j in a:
          Fe_xiv = ch.ion('fe_14',temperature=T_sec3[:,i,j],eDensity=n_sec3[:,i,j],abundance='sun_coronal_1992_feldman_ext')  
          Fe_xiv.intensity()
          wavelengths=Fe_xiv.Intensity['wvl']
          intensities=Fe_xiv.Intensity['intensity'] 
          index_ord=np.argsort(wavelengths)
          wvl_ord=wavelengths[index_ord]
          intensities_ord=intensities[:,index_ord]

          i_211=find_closest(wvl_ord,211)
          r_wvl=np.arange(i_211-5,i_211+6)
          r_int=intensities_ord[:,r_wvl]
            
          #Some conditions
          g_211=np.sum(r_int,axis=1)
          g_211[T_sec3[:,i,j]<1.e5]=0
          g_211[n_sec3[:,i,j]>1.e11]=0

          #Integral with the trapezoidal rule for a non-uniform spacing grid
          n_H=n_sec3[:,i,j]/1.2
          integrando=g_211*n_H*n_sec3[:,i,j]
          integral_i=np.trapz(integrando, x=z_sec)
          integral_map.append(integral_i)

  integral_map=np.reshape(integral_map,(resolution,resolution))

  #Saving the fits
  hdu = fits.PrimaryHDU(integral_map)
  hdu.writeto('integral_map_4.fits',overwrite=True)
