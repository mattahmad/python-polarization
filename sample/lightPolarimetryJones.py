import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix, create_Jones_matrices
from py_pol.stokes import Stokes, create_Stokes
from py_pol.mueller import Mueller, create_Mueller

def light_polarimeter_measurement(E, angleR, angleD):
    """This function simulates the light polarimeter measurements.

    Parameters:
        E: Incident light wave Jones vector.
        angleR: Rotation angle of the quarter-wave plate.
        angleD: Rotating angle of the diattenuator (polarizer).

    Returns:
        I: Measured intensity."""
    # Rotate the optical elements
    Jr_rotated = Jr.rotate(angle=angleR, keep=True)
    Jd_rotated = Jd.rotate(angle=angleD, keep=True)
    # Multiply all objects
    E_final = Jd_rotated * Jr_rotated * E
    # Measure the intensity
    I = E_final.parameters.intensity()
    # Return
    return I

# Start by defining the light wave
I_wave = np.random.rand() * 4 + 1
alpha_wave = np.random.rand() * 90*degrees
delay_wave = np.random.rand() * 360*degrees
E_wave = Jones_vector('Light wave')
E_wave.general_charac_angles(intensity=I_wave, alpha=alpha_wave, delay=delay_wave)
print(E_wave)

# Now, define the optical elements
Jr, Jd = create_Jones_matrices(('Quarter-wave plate', 'Perfect polarizer'))
Jr.quarter_waveplate(azimuth=0)
Jd.diattenuator_perfect(azimuth=0)

# Now, define the rotation angles for the optical elements
angles_R = np.array([0*degrees, 0*degrees, 135*degrees, 45*degrees])
angles_D = np.array([0*degrees, 90*degrees, 0*degrees, 45*degrees])

# Make the measurements
If1, If2, If3, If4 = light_polarimeter_measurement(E_wave, angles_R, angles_D)
print(If1, If2, If3, If4)

# Calculate the parameters
I_calc = If1 + If2
alpha_calc = np.arctan(np.sqrt(If2/If1))
delay_calc = np.arctan2(2*If3 - I_calc, 2*If4 - I_calc) % (2*np.pi)

# Compare the result
print('Comparison')
print('  - Intensity:')
print('      o   Original:   ', I_wave)
print('      o   Measured:   ', I_calc)
print('  - Alpha (deg):')
print('      o   Original:   ', alpha_wave/degrees)
print('      o   Measured:   ', alpha_calc/degrees)
print('  - Delay (deg):')
print('      o   Original:   ', delay_wave/degrees)
print('      o   Measured:   ', delay_calc/degrees)

def intensity_difference(parameters):
    """This function calculates the difference between the measurements and the model.

    Parameters:
        parameters: List containing the intensity, the alpha and the delay.

    Returns:
        dI: Intensity difference."""
    # Measure the intensity
    I_meas = light_polarimeter_measurement(E_wave, angles_R, angles_D)
    # Calculate the intensity from the model
    E_model = Jones_vector('Test wave')
    E_model.general_charac_angles(intensity=parameters[0], alpha=parameters[1], delay=parameters[2])
    I_model = light_polarimeter_measurement(E_model, angles_R, angles_D)
    # Calculate the difference and return it
    dI = I_meas - I_model
    return dI

def polarimetry_fit(fun, N, tol=None, verbose=True):
    """Function that calculates the parameters of the Jones vector using N measurements and a fitting algorithm.

    Parameters:
        N (int): Number of measurements.
        tol (float): Tolerance of the fitting algorithm.
        verbose (bool): If True, the function prints the comparison between real and fit variables.

    Returns.
        x (list): List containing the fit parameters."""
    # Prepare the fitting
    sup_limit = np.array([5, 90*degrees, 360*degrees])
    inf_limit = np.zeros(3)
    x0 = np.random.rand(3) * sup_limit
    # Make the fitting
    result = least_squares(fun=fun, x0=x0, bounds=(inf_limit, sup_limit), ftol=tol, xtol=tol, gtol=tol)

    # Compare the result
    if verbose:
        print('Comparison')
        print('  - Intensity:')
        print('      o   Original:   ', I_wave)
        print('      o   Measured:   ', result.x[0])
        print('  - Alpha (deg):')
        print('      o   Original:   ', alpha_wave/degrees)
        print('      o   Measured:   ', result.x[1]/degrees)
        print('  - Delay (deg):')
        print('      o   Original:   ', delay_wave/degrees)
        print('      o   Measured:   ', result.x[2]/degrees)

    # Return
    return result.x

# Start by defining the light wave
I_wave = np.random.rand() * 5
alpha_wave = np.random.rand() * 90*degrees
delay_wave = np.random.rand() * 360*degrees
E_wave = Jones_vector('Light wave')
E_wave.general_charac_angles(intensity=I_wave, alpha=alpha_wave, delay=delay_wave)

# Make the experiment
N = 20
tol = 1e-12
polarimetry_fit(fun=intensity_difference, N=N, tol=tol);