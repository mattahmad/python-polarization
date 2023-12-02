import numpy as np
import matplotlib.pyplot as plt

from py_pol import degrees
from py_pol.stokes import Stokes
from py_pol.mueller import Mueller

S0 = Stokes('Incident wave')
S0.linear_light(intensity=1, azimuth=0*degrees)

Nmax = 30
Narray = np.arange(2,Nmax+1)
J = Mueller('System')
I = np.zeros(Nmax-1)

# For loop to calculate the intensities
for ind, N in enumerate(Narray):
    # Create the system of polarizers
    angles = np.linspace(90*degrees, 0*degrees, N) # Start by 90º and end in 0º because the rightest  element is the first to be crossed by the light wave
    J.diattenuator_linear(Tmax=0.99, Tmin=1e-4, azimuth=angles)
    # Multiply them
    J.prod()
    # Multiply the system matrix by the incident light wave to obtain the output
    S_final = J * S0
    # Store the intensity
    I[ind] = S_final.parameters.intensity(out_number=False)

    # Plot the result
plt.figure()
plt.plot(Narray, I, 'ko')
plt.xlabel('# polarizers', fontsize=22)
plt.ylabel('$I_{transmited}$', fontsize=22)
plt.title('Imperfect polarizers', fontsize=24)

# Find the maximum
Imax = np.max(I)
Nmax = np.argmax(I) + 2
print('The maximum obtained intensity is {} for {} polarizers.'.format(Imax, Nmax))