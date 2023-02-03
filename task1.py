#Add this part to make it work on windows
import os
os.environ['XUVTOP']='C:\\Users\\sarit\\CHIANTI'
os.environ["HOME"]="C:\\Users\\sarit\\CHIANTI"

#Packages
import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt


#Defining temperature and electron density
temp=np.logspace(4.69897,7.69897,46) #kelvin
edens=3.e+9 #cm-3

#Gain function
Fe_xviii = ch.ion('fe_18',temperature=temp,eDensity=edens,abundance='sun_coronal_1992_feldman_ext')  

#Level transitions within the Fe XVIII ion that lead to line emission
Fe_xviii.intensity() 

'''
#Intensity is a dictionary with the following keys
for ii in Fe_xviii.Intensity.keys(): print(ii)
#To see the size and shape of the values: 
for ii in Fe_xviii.Intensity.values(): print(ii.shape)
'''

wavelengths=Fe_xviii.Intensity['wvl'] #the wavelengths are not ordered
#Gain functions
intensities=Fe_xviii.Intensity['intensity']#[Temperature,wavelengths], gain for different temperatures and wavelengths
#Sorting the wavelength in ascending order and the gain
index_ord=np.argsort(wavelengths)
wvl_ord=wavelengths[index_ord] 
intensities_ord=intensities[:,index_ord]

#Finding maximums
ind_ord = np.unravel_index(np.argmax(intensities_ord, axis=None), intensities_ord.shape) 
index_Tmax=ind_ord[0]
index_wvlmax=ind_ord[1]


#---Task 1a---
#Plot the intensities as a function of wavelength
plt.figure(1)
markerline, stemlines, baseline = plt.stem(wvl_ord, intensities_ord[index_Tmax,:], linefmt='k', markerfmt='D')
markerline.set_visible(False)
baseline.set_visible(False)
plt.xlim(90,96)
plt.ylim(0,1.40e-25)
plt.xlabel(r'Wavelength ($\AA$)')
plt.ylabel(r'Gain (erg $cm^{-2} s^{-1} sr^{-1}$)')
plt.title(r'Fe XVII T=7.9e06(K) Density=3.0e+09($cm^{-3}$)')
plt.text(89.017-0.3,0.03e-25, 'Fe XVIII 89.017', c = 'k', fontsize=9,rotation=90)
plt.text(89.722-0.3,0.03e-25, 'Fe XVIII 89.722', c = 'k', fontsize=9,rotation=90)
plt.text(93.932-0.8,0.03e-25, 'Fe XVIII 93.932', c = 'k', fontsize=9,rotation=90)
plt.text(99.565-0.3,0.03e-25, 'Fe XVIII 99.565', c = 'k', fontsize=9,rotation=90)
plt.text(100.735-0.3,0.03e-25, 'Fe XVIII 100.735', c = 'k', fontsize=9,rotation=90)
plt.text(103.948-0.8,0.03e-25, 'Fe XVIII 103.948', c = 'k', fontsize=9,rotation=90)
plt.text(104.324-0.3,0.03e-25, 'Fe XVIII 104.324', c = 'k', fontsize=9,rotation=90)
plt.text(105.475-0.3,0.03e-25, 'Fe XVIII 105.475', c = 'k', fontsize=9,rotation=90)
plt.text(105.987-0.3,0.03e-25, 'Fe XVIII 105.987', c = 'k', fontsize=9,rotation=90)
plt.text(107.235-0.3,0.03e-25, 'Fe XVIII 107.235', c = 'k', fontsize=9,rotation=90)
plt.text(109.606-0.3,0.03e-25, 'Fe XVIII 109.606', c = 'k', fontsize=9,rotation=90)
plt.show()
#plt.savefig('task1a.eps',format='eps')




#---Task 1b---

gain_total=intensities_ord.sum(axis=1)

plt.figure(3)

plt.loglog(temp,gain_total)
plt.xlim(6e+5,1e+8)
plt.xlabel(r'Temperature (K)')
plt.ylabel(r'Gain$_{total}$ (erg $cm^{-2} s^{-1} sr^{-1}$)')


#---Task 1c---
plt.figure(4)

edens=np.array([3e+7,3e+8,3e+9,3e+10])
for i in edens:
    Fe_xviii = ch.ion('fe_18',temperature=temp,eDensity=i,abundance='sun_coronal_1992_feldman_ext')  
    Fe_xviii.intensity()
    wavelengths=Fe_xviii.Intensity['wvl'] 
    intensities=Fe_xviii.Intensity['intensity'] ##gain functions, [Temperature,wavelengths] sets de ganancia para diferentes temperaturas y longitudes de onda
    index_ord=np.argsort(wavelengths)
    wvl_ord=wavelengths[index_ord] # ordenamos de forma ascendente
    intensities_ord=intensities[:,index_ord] 
    gain_total=intensities_ord.sum(axis=1)
    plt.loglog(temp,gain_total,label=r"$n_e$={:.2e}".format(i))


plt.loglog(temp,gain_total)
plt.xlim(6e+5,1e+8)
plt.title(r'g$_{total}$(T,$n_e$) for FeXVIII')
plt.xlabel(r'Temperature (K)')
plt.ylabel(r'Gain$_{total}$ (erg $cm^{-2} s^{-1} sr^{-1}$)')
plt.legend()

#plt.savefig('task1c.eps',format='eps')
