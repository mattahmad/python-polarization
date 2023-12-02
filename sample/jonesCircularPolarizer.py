import numpy as np

from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix

# Incident wave
E0 = Jones_vector('Incident wave')
E0.circular_light(intensity=1)
print(E0)