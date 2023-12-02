import numpy as np

from py_pol import degrees
from py_pol.stokes import Stokes
from py_pol.mueller import Mueller

# Incident wave
S0 = Stokes('Incident wave')
S0.circular_light(intensity=1)
print(S0)