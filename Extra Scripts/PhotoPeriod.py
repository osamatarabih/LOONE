# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 14:07:18 2023

@author: osama
"""

# Import modules
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# Define function
def photoperiod(phi,doy,verbose=False):
    
    phi = np.radians(phi) # Convert to radians
    light_intensity = 2.206 * 10**-3

    C = np.sin(np.radians(23.44)) # sin of the obliquity of 23.44 degrees.
    B = -4.76 - 1.03 * np.log(light_intensity) # Eq. [5]. Angle of the sun below the horizon. Civil twilight is -4.76 degrees.

    # Calculations
    alpha = np.radians(90 + B) # Eq. [6]. Value at sunrise and sunset.
    M = 0.9856*doy - 3.251 # Eq. [4].
    lmd = M + 1.916*np.sin(np.radians(M)) + 0.020*np.sin(np.radians(2*M)) + 282.565 # Eq. [3]. Lambda
    delta = np.arcsin(C*np.sin(np.radians(lmd))) # Eq. [2].

    # Defining sec(x) = 1/cos(x)
    P = 2/15 * np.degrees( np.arccos( np.cos(alpha) * (1/np.cos(phi)) * (1/np.cos(delta)) - np.tan(phi) * np.tan(delta) ) ) # Eq. [1].

    # Print results in order for each computation to match example in paper
    if verbose:
        print('Input latitude =', np.degrees(phi))
        print('[Eq 5] B =', B)
        print('[Eq 6] alpha =', np.degrees(alpha))
        print('[Eq 4] M =', M[0])
        print('[Eq 3] Lambda =', lmd[0])
        print('[Eq 2] delta=', np.degrees(delta[0]))
        print('[Eq 1] Daylength =', P[0])
    
    return P
# Invoke function with scalars
# phi = 26.982052;  # Latitude for consistency with notation in literature.
# doy = np.array([201]); # Day of the year. Julian calendar. Day from January 1.
# P = photoperiod(phi,doy,verbose=True)
# print('Photoperiod: ' + str(np.round(P[0],2)) + ' hours/day')
# Multiple inputs call
phi = 26.982052;
doy = np.arange(1,365);
P = photoperiod(phi,doy)
plt.figure(figsize=(8,6))
plt.plot(doy,P)
plt.title('Latitude:' + str(phi))
plt.xlabel('Day of the year', size=14)
plt.ylabel('Photoperiod (hours per day)', size=14)
plt.show()
print(P)
PhotoPeriod = pd.DataFrame()
PhotoPeriod['Day'] = doy
PhotoPeriod['PhotoPeriod'] = P
PhotoPeriod.to_csv('C:/Work/Research/LOONE/Nitrogen Module/PhotoPeriod.csv')
