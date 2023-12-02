import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix, create_Jones_matrices
from py_pol.stokes import Stokes, create_Stokes
from py_pol.mueller import Mueller, create_Mueller

def Mueller_polarimeter_measurement(Msample, Mpol, S_ini, angles, Ierror=0):
    """This function simulates the light polarimeter measurements.

    Parameters:
        Msample (Mueller): Mueller matrix of the sample.
        Mpol (list): List with the four Mueller objects of the polarimeter: Md1, Mr1, Mr2, Md2.
        Sini (Stokes): Stokes object of the illumination.
        angles (numpy.ndarray): 4xN array with the angles used to perform the measurement.
        Ierror (float): Intensity error amplitude.

    Returns:
        Ik: Measured intensities vector.
        b: b matrix."""
    # Extract the information from the input parameters
    N = angles.shape[1]
    Md1, Mr1, Mr2, Md2 = Mpol
    # Rotate the optical elements
    Md1_rotated = Md1.rotate(angle=angles[0,:], keep=True)
    Mr1_rotated = Mr1.rotate(angle=angles[1,:], keep=True)
    Mr2_rotated = Mr2.rotate(angle=angles[2,:], keep=True)
    Md2_rotated = Md2.rotate(angle=angles[3,:], keep=True)
    # Multiply the objects as required
    Sg = Mr1_rotated * Md1_rotated * S_ini
    A = Md2_rotated * Mr2_rotated
    S_final = A * Msample * Sg
    # Measure the intensity
    Ik = S_final.parameters.intensity() + np.random.randn(N) * Ierror
    # Calculate the b matrix
    Ak = A.parameters.matrix()   # Extract the 4x4xN matrix of the Mueller object
    Ak = Ak[0,:,:]   # Extract all the first rows
#     Ak = np.transpose(Ak)   # Transpose so Ak rows corrrespond to different measurements
    Sk = Sg.parameters.matrix()   # Extract the 4xN matrix of the Stokes object
#     Sk = np.transpose(Sk)   # Transpose so Sk rows corrrespond to different measurements
    b = np.zeros((N,16))
    for ind in range(N):
        b[ind,:] = np.outer(Ak[:,ind], Sk[:,ind]).flatten()
    # Return
    return Ik, b

# Start by defining the polarimeter system
S_ini = Stokes('Light source')
S_ini.circular_light(intensity=1)
Mpol = create_Mueller(('Diattenuator 1', 'Retarder 1', 'Retarder 2', 'Diattenuator 2'))
Mpol[0].diattenuator_perfect(azimuth=0)
Mpol[1].quarter_waveplate(azimuth=0)
Mpol[2].quarter_waveplate(azimuth=0)
Mpol[3].diattenuator_perfect(azimuth=0)

# Create the sample parameters
p1 = np.random.rand() * 0.5 + 0.5
p2 = np.random.rand() * 0.5
alphaD = np.random.rand() * 90*degrees
delayD = np.random.rand() * 360*degrees

R = np.random.rand() * 180*degrees
alphaR = np.random.rand() * 90*degrees
delayR = np.random.rand() * 360*degrees

d = [np.random.rand(), np.random.rand(), np.random.rand()]

# Create the sample
Md, Mr, Mp = create_Mueller(N=3)
Md.diattenuator_charac_angles(p1=p1, p2=p2, alpha=alphaD, delay=delayD)
Mr.retarder_charac_angles(R=R, alpha=alphaR, delay=delayR)
Mp.depolarizer_diagonal(d=d)

Msample = Md*Mr*Mp
Msample.name = 'Sample'

# Now, define the rotation angles for the optical elements
angles = np.random.rand(4,16) * 180*degrees

# Make the measurements
Ik, b = Mueller_polarimeter_measurement(Msample=Msample, Mpol=Mpol, S_ini=S_ini, angles=angles)

# Calculate the original Mueller matrix
b_inv = np.linalg.inv(b)
Ik = np.reshape(Ik, (16,1))
M_calc = Mueller('Calculated')
M_calc.from_matrix(M=b_inv@Ik)

# Compare the result
print(Msample, M_calc)