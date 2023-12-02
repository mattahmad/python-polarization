import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt

from py_pol import degrees
from py_pol.jones_vector import Jones_vector
from py_pol.jones_matrix import Jones_matrix, create_Jones_matrices
from py_pol.stokes import Stokes, create_Stokes
from py_pol.mueller import Mueller, create_Mueller

def light_polarimeter_measurement_Mueller(S, angleR, angleD, Ierror=0):
    """This function simulates the light polarimeter measurements.

    Parameters:
        S: Incident light wave Stokes vector.
        angleR: Rotation angle of the quarter-wave plate.
        angleD: Rotating angle of the diattenuator (polarizer).

    Returns:
        Ik: Measured intensities vector.
        Ak: Matrix with the first rows of the analizer matrices."""
    N = angleR.size
    # Rotate the optical elements
    Mr_rotated = Mr.rotate(angle=angleR, keep=True)
    Md_rotated = Md.rotate(angle=angleD, keep=True)
    A = Md_rotated * Mr_rotated
    # Multiply all objects
    S_final = A * S
    # Measure the intensity
    Ik = S_final.parameters.intensity() + np.random.randn(N) * Ierror
    # Calculate the Ak matrix
    Ak = A.parameters.matrix()   # Extract the 4x4xN matrix of the py_pol object
    Ak = Ak[0,:,:]   # Extract all the first rows
    Ak = np.transpose(Ak)   # Transpose so each row of the final matrix corresponds to a first row of the analizer matrices
    # Return
    return Ik, Ak

# Start by defining the light wave
I_wave = np.random.rand(1) * 5
alpha_wave = np.random.rand(1) * 90*degrees
delay_wave = np.random.rand(1) * 360*degrees
degree_pol = np.random.rand(1)
S_wave = Stokes('Original')
S_wave.general_charac_angles(intensity=I_wave, alpha=alpha_wave, delay=delay_wave, degree_pol=degree_pol)

# Now, define the optical elements
Mr, Md = create_Mueller(('Quarter-wave plate', 'Perfect polarizer'))
Mr.quarter_waveplate(azimuth=0)
Md.diattenuator_perfect(azimuth=0)

# Now, define the rotation angles for the optical elements
angles_R = np.random.rand(4) * 180*degrees
angles_D = np.random.rand(4) * 180*degrees

# Make the measurements
Ik, Ak = light_polarimeter_measurement_Mueller(S_wave, angles_R, angles_D)

# Calculate the original Stokes vector
Ak_inv = np.linalg.inv(Ak)
Ik = np.reshape(Ik, (4,1))
S_calc = Stokes('Calculated')
S_calc.from_matrix(M=Ak_inv@Ik)

# Compare the result
print(S_wave, S_calc)

def light_polariemter_experiment(N, Ierror):
    """This funcion simulates a polarimetry experiment with errors in the detection.

    Parameters:
        N (int): Number of measurements.
        Ierror (float): Intensity error amplitude.

    Returns:
        dif (numpy.ndarray): Difference between the calculated and the original array."""
    # Start by defining the light wave
    alpha_wave = np.random.rand(1) * 90*degrees
    delay_wave = np.random.rand(1) * 360*degrees
    degree_pol = np.random.rand(1)
    S_wave = Stokes('Original')
    S_wave.general_charac_angles(intensity=I_wave, alpha=alpha_wave, delay=delay_wave, degree_pol=degree_pol)

    # Now, define the rotation angles for the optical elements
    angles_R = np.random.rand(N) * 180*degrees
    angles_D = np.random.rand(N) * 180*degrees

    # Make the measurements
    Ik, Ak = light_polarimeter_measurement_Mueller(S_wave, angles_R, angles_D, Ierror=Ierror)

    # Calculate the Stokes vector
    Ak_T = np.transpose(Ak)
    Ak_inv = np.linalg.inv(Ak_T @ Ak)
    Ak_inv = Ak_inv @ Ak_T
    Ik = np.reshape(Ik, (N,1))
    S_calc = Ak_inv @ Ik

    # Calculate the difference
    dif = np.squeeze(S_wave.M - S_calc)

    # Return
    return dif

# Define some variables
Nstat = 100      # Number of experiments done under the same condition.
factors = np.array([0, 0.005, 0.01, 0.02])    # Array of error factors
Ns = np.arange(10, 100, 10)    # Array of the number of measurements
tol = 1e-12
I_wave = 1    # We will maintain this parameter constant so we can compare the errors easily

# Create the result variables
error = np.zeros((4, Ns.size, factors.size))
legend = []

# Make the loops changing the variables
for indI, factor in enumerate(factors):   # Loop in the error amplitude
    for indN, N in enumerate(Ns):   # Loop in the number of measurements
        # Create the temporal variables
        aux = np.zeros((4, Nstat))
        # Repeat the same conditions several times
        for indS in range(Nstat):
            # Do the experiment
            dif = light_polariemter_experiment(N=N, Ierror=I_wave*factor);
            # Calculate the errors
            aux[:,indS] = dif
        # Calculate the mean square error
        error[:, indN, indI] = np.linalg.norm(aux, axis=1) / Nstat
    # Create the legend for the plots
    legend += ['Ierror = {} %'.format(factor*100)]

    # Plot the results
plt.figure(figsize=(12,12))
plt.subplot(2,2,1)
plt.plot(Ns, error[0,:,:])
plt.legend(legend)
plt.xlabel('Number of measurements')
plt.ylabel('S0')

plt.subplot(2,2,2)
plt.plot(Ns, error[1,:,:])
plt.legend(legend)
plt.xlabel('Number of measurements')
plt.ylabel('S1')

plt.subplot(2,2,3)
plt.plot(Ns, error[2,:,:])
plt.legend(legend)
plt.xlabel('Number of measurements')
plt.ylabel('S2')

plt.subplot(2,2,4)
plt.plot(Ns, error[3,:,:])
plt.legend(legend)
plt.xlabel('Number of measurements')
plt.ylabel('S3')

plt.show()

