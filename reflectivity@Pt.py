#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 14:34:04 2023

@author: bruno.souza
"""

####### Calculate the reflectivity of the material, as well as the total power and the filtered power #######
# this script obtain the graphics:
    # Reflectivity vs Energy
    # Total flux vs Energy 
    # Filtered flux vs Energy 


# When the beam passes through the mirror, part of its power is absorbed. 
# Here, we are referring to the absorbed portion as filtered power.
    
# Packages and functions
from optlnls.mirror import reflectivity_xrays
from scipy import integrate as inte
from scipy import constants as con
import numpy as np
import matplotlib.pyplot as plt


# Parameters
mat  = 'Pt'
dens = 21.45 # material density [g/cm3]
atm  = 195.08 # atomic mass
angs = 2.25e-3 # surface angle    [rad]
angn = 90-np.rad2deg(angs)
Eo   = 100   #  [eV] 
Ef   = 80e3 #   [eV]
n    = 10000 # [dimensionless]
e_charge=con.elementary_charge  # elementary charge [C]
e_vector=(np.linspace(Eo,Ef,n)) # energy vector     [eV]


# Reflectivity calculate
refl_vector= reflectivity_xrays(material=mat, density=dens, atomic_mass=atm, energy_eV=e_vector, angle_normal_deg=angn)
trans_vector= 1-refl_vector

# Create a vector with data
vetor_data= np.array([e_vector,refl_vector]).transpose()

# Save .txt
filename = 'Refl_Pt_E100_400e3.txt'
np.savetxt(filename, vetor_data, fmt='%.6f', delimiter='\t', header="Energy (eV)\tReflectivity")

# Plot Reflectivity Curve of Pt
plt.figure()
plt.title("Pt Reflectivity")
plt.xlabel("Energy [keV]")
plt.ylabel("Reflectivity")
plt.xlim(0,Ef/1000)
plt.ylim(0,1)
plt.grid()
plt.plot(e_vector/1000,refl_vector, ':')
# Save png 
plt.savefig(str(mat)+'_Reflectivity', dpi=600)

# Open .txt from Spectra -> create 3 vectors
e_txt,flu_txt,fil_txt=np.genfromtxt('Flux_QUATI_M1_E1.txt',unpack=True,usecols=(0,1,5),skip_header=2)
# e_txt   [eV]
# flu_txt [ph/sec/0.1%BW]
# fil_txt [ph/sec/0.1%BW]

### Plot Flux M1
plt.figure()
plt.title("Flux M1 - QUATI")
plt.xlabel("Energy [keV]")
plt.ylabel("FLux [ph/sec/0.1%BW/100mA]")
plt.grid(True, which="both", ls="-")
plt.plot(e_txt/1000,flu_txt, '--',label='Total Flux')
plt.plot(e_txt/1000,fil_txt, '-',label='Filtered by M1')
plt.legend()
plt.xscale('log')
plt.yscale('log')

# Save png 
plt.savefig('Flux.png', dpi=600)


### Unit conversion [ph/sec/0.1%BW] -> [ph/sec/eV]
flux_ev=flu_txt*(1000/e_txt)
fil_ev=fil_txt*(1000/e_txt)

### Taking into account the energy of each photon
a=flux_ev*(e_txt)*(e_charge)   #  [J.ph/sec/ev]
b=fil_ev*(e_txt)*(e_charge)    #  [J.ph/sec/ev]

### Integrate
total_power=inte.simps(a,e_txt)     # [W]
filtered_power=inte.simps(b,e_txt)  # [W]

print('\n')
print('Total Power = %.3f W' %total_power)
print('Filtered Power = %.3f W' %filtered_power)

total_flux=inte.simps(flux_ev,e_txt)
filtered_flux=inte.simps(fil_ev,e_txt)

print('\n')
print('Total Flux = %.3e ph/sec/100mA' %total_flux)
print('Filtered Flux = %.3e ph/sec/100mA' %filtered_flux)


