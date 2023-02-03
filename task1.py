

#Agregar esta parte para que funcione en windows
import os

os.environ['XUVTOP']='C:\\Users\\sarit\\CHIANTI'
os.environ["HOME"]="C:\\Users\\sarit\\CHIANTI"



import ChiantiPy.core as ch
import numpy as np
import matplotlib.pyplot as plt


plt.close('all')
temp=np.logspace(4.69897,7.69897,46) #kelvin
#temp=temp[33]
#edens=np.ones(46)*3.e+9
edens=3.e+9 #cm-3

Fe_xviii = ch.ion('fe_18',temperature=temp,eDensity=edens,abundance='sun_coronal_1992_feldman_ext')  

#Fe_xviii.popPlot()
#Fe_xviii.intensityPlot(wvlRange=[210.,220])
#plt.show()

Fe_xviii.intensity() #level transitions within the Fe XVIII ion that lead to line emission


wavelengths=Fe_xviii.Intensity['wvl'] #las longitudes de onda no estan ordenadas 
#da como salida una matriz con las temperaturas que hemos definido y las longitudes de onda donde se han encontrado los niveles de transición
#para las el ion
intensities=Fe_xviii.Intensity['intensity'] ##gain functions, [Temperature,wavelengths] sets de ganancia para diferentes temperaturas y longitudes de onda
index_ord=np.argsort(wavelengths)
wvl_ord=wavelengths[index_ord] # ordenamos de forma ascendente
intensities_ord=intensities[:,index_ord] #ordenamos la ganacia en longitudes de onda


ind_ord = np.unravel_index(np.argmax(intensities_ord, axis=None), intensities_ord.shape) #posicion del valor maximo
index_Tmax=ind_ord[0]
index_wvlmax=ind_ord[1]


#Task 1a
plt.figure(1)
#plt.stem(wvl_ord , intensities_ord[index_Tmax,:],linefmt='k',markerfmt='D',markersize=0.1)
markerline, stemlines, baseline = plt.stem(wvl_ord, intensities_ord[index_Tmax,:], linefmt='k', markerfmt='D')
#markerline.set_color('k')
#markerline.set_markersize(2)
markerline.set_visible(False)
baseline.set_visible(False)
#plt.plot(wvl_ord , intensities_ord[index_Tmax,:])
plt.xlim(90,96)
plt.ylim(0,1.40e-25)
#plt.yscale("log") 
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

#plt.figure(2)

#Fe_xviii.intensityPlot(wvlRange=[85.,112.]) #esto solo printea las 10 más intensas y una T random

'''
#Intensity is a dictionary with the following keys
for ii in Fe_xviii.Intensity.keys(): print(ii)

#To see the size and shape of the values: 
for ii in Fe_xviii.Intensity.values(): print(ii.shape)
'''


#Task 1b, Task 1c

gain_total=intensities_ord.sum(axis=1)

plt.figure(3)

plt.loglog(temp,gain_total)
plt.xlim(6e+5,1e+8)
plt.xlabel(r'Temperature (K)')
plt.ylabel(r'Gain$_{total}$ (erg $cm^{-2} s^{-1} sr^{-1}$)')


#Task 1c
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
