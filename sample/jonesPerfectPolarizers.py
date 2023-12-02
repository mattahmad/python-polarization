import numpy as np

from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix

# First polarizer
P0 = Jones_matrix('Polarizer 0')
P0.diattenuator_perfect(azimuth=50 * degrees)
print(P0)

# Second polarizer
P1 = Jones_matrix('Polarizer 1')
angles = np.linspace(0, 360*degrees, 361) # Steps of 1 degree
P1.diattenuator_perfect(azimuth=angles);
print(P1[0], P1[180], P1[-1])