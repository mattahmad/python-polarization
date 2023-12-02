import numpy as np

from py_pol import degrees
from py_pol.stokes import Stokes
from py_pol.mueller import Mueller

# Incident wave
S0 = Stokes('Incident wave')
S0.circular_light(intensity=1)

# First polarizer
P0 = Mueller('Polarizer 0')
P0.diattenuator_perfect(azimuth=50 * degrees)

# Second polarizer
P1 = Mueller('Polarizer 0')
angles = np.linspace(0, 360*degrees, 361) # Steps of 1 degree
P1.diattenuator_perfect(azimuth=angles);

S_final = P1 * P0 * S0
S_final.name = 'Output wave'

I_perfect = S_final.parameters.intensity(draw=True)