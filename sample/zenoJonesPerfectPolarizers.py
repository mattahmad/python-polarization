from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix

import matplotlib.pyplot as plt
import numpy as np

E0 = Jones_vector('Incident wave')
E0.linear_light(intensity=1, azimuth=0*degrees)

Nmax = 30
Narray = np.arange(2,Nmax+1)
J = Jones_matrix('System')
I = np.zeros(Nmax-1)

# For loop to calculate the intensities
for ind, N in enumerate(Narray):
    # Create the system of polarizers
    angles = np.linspace(90*degrees, 0*degrees, N) # Start by 90º and end in 0º because the rightest  element is the first to be crossed by the light wave
    J.diattenuator_perfect(azimuth=angles)
    # Multiply them
    J.prod()
    # Multiply the system matrix by the incident light wave to obtain the output
    E_final = J * E0
    # Store the intensity
    I[ind] = E_final.parameters.intensity(out_number=False)

    # Plot the result
plt.figure()
plt.plot(Narray, I, 'ko')
plt.xlabel('# polarizers', fontsize=22)
plt.ylabel('$I_{transmited}$', fontsize=22)
plt.title('Perfect polarizers', fontsize=24);