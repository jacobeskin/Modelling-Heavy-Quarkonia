#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

# Read the wavefunction from a text file into a NumPy array
filename = raw_input("Input file name: ")
wf = np.loadtxt(filename)

# Create the radial wavefunction from the reduced wavefunction
R = np.zeros(wf.size)
for i in range(0, wf.size):
    R[i] = wf[i]/(i+1)

# Plot the reduce radial wavefunctio and the radial wavefunction
x = np.arange(0, wf.size)

plt.figure(1)

plt.subplot(211)
plt.title('Reduced radial wavefunction u(r)', fontweight='bold', fontsize=17)
plt.xlabel('r (am)', fontweight='bold', fontsize=15)
plt.ylabel('u(r)', fontweight='bold', fontsize=15)
plt.plot(x, wf, 'b', linewidth=2.0)
plt.axhline(0, color='black')

plt.subplot(212)
plt.title('Radial wavefunction R(r)', fontweight='bold', fontsize=17)
plt.xlabel('r (am)', fontweight='bold', fontsize=15)
plt.ylabel('R(r)', fontweight='bold', fontsize=15)
plt.plot(x, R, 'r', linewidth=2.0)
plt.axhline(0, color='black')

plt.show()




