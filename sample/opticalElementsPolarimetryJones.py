import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix, create_Jones_matrices
from py_pol.stokes import Stokes, create_Stokes
from py_pol.mueller import Mueller, create_Mueller

def sample_polarimeter_measurement(Jsample, angleDG, angleRG, angleRA, angleDA):
    """This function simulates the light polarimeter measurements.

    Parameters:
        E: Incident light wave Jones vector.
        angleDG: Rotating angle of the diattenuator (polarizer) of the PSG.
        angleRG: Rotation angle of the quarter-wave plate of the PSG.
        angleRA: Rotation angle of the quarter-wave plate of the PSA.
        angleDA: Rotating angle of the diattenuator (polarizer) of the PSA.

    Returns:
        I: Measured intensity."""
    # Rotate the optical elements
    E = Jones_vector().circular_light(intensity=2)   # Intensity = 2 produces a normalized light state of intensity 1 after the PSG
    Jdg = Jones_matrix().diattenuator_perfect(azimuth=angleDG)
    Jrg = Jones_matrix().quarter_waveplate(azimuth=angleRG)
    Jra = Jones_matrix().quarter_waveplate(azimuth=angleRA)
    Jda = Jones_matrix().diattenuator_perfect(azimuth=angleDA)
    # Multiply all objects
    E_final = Jda * Jra * Jsample * Jrg * Jdg * E
    # Measure the intensity
    I = E_final.parameters.intensity()
    # Return
    return I

def matrix_parameters(I):
    """This function calculates the parameters of the Jones matrix from the intensity measurements.

    Parameters:
        I: Array of intensities.

    Returns:
        J0, J1, J2, J3, d1, d2, d3: Matrix parameters."""
    # Calculate modules
    J0, J1, J2, J3 = np.sqrt(I[:4])
    # Calculate phases
    d1 = np.arctan2(2*I[4] - I[1] - I[0], 2*I[5] - I[1] - I[0]) % (360*degrees)
    d2 = np.arctan2(2*I[6] - I[2] - I[0], 2*I[7] - I[2] - I[0]) % (360*degrees)
    d3 = (np.arctan2(2*I[8] - I[2] - I[3], 2*I[9] - I[2] - I[3]) + d2) % (360*degrees)
    # Return
    return J0, J1, J2, J3, d1, d2, d3

# Generate a random Jones matrix for the sample
J0, J1, J2, J3 = np.random.rand(4)
d1, d2, d3 = np.random.rand(3) * 360*degrees
Jsample = Jones_matrix().from_components([J0, J1*np.exp(1j*d1), J2*np.exp(1j*d2), J3*np.exp(1j*d3)])

# Calculate the intensities
angleDG = np.array([0, 90, 0, 90, 45, 45, 0, 0, 45, 45]) * degrees
angleRG = np.array([0, 90, 0, 90, 0, 45, 0, 0, 0, 45]) * degrees
angleRA = np.array([0, 0, 90, 90, 0, 0, 0, 45, 90, 90]) * degrees
angleDA = np.array([0, 0, 90, 90, 0, 0, 45, 45, 90, 90]) * degrees
I = sample_polarimeter_measurement(Jsample, angleDG, angleRG, angleRA, angleDA)

# Calculate the matrix parameters
J0_calc, J1_calc, J2_calc, J3_calc, d1_calc, d2_calc, d3_calc = matrix_parameters(I)

# Compare results
print("Comparison between original and calculated parameters:")
original = [J0, J1, J2, J3, d1/degrees, d2/degrees, d3/degrees]
calculated = [J0_calc, J1_calc, J2_calc, J3_calc, d1_calc/degrees, d2_calc/degrees, d3_calc/degrees]
parameter = ["J0", "J1", "J2", "J3", "d1 (deg)", "d2 (deg)", "d3 (deg)"]
for ind in range(7):
    print("- {}:   Original = {:.4f};    Calculated = {:.4f}".format(parameter[ind], original[ind], calculated[ind]))
