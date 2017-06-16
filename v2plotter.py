import matplotlib.pyplot as plt
import numpy as np

# Creates a scatterplot of all the <v^2> values for Charmonium and Bottomonium

# The data

CS = [0.2974, 0.2917, 0.3504, 0.3706, 0.4404, 0.4303, 0.5134, 0.5064]
CP = [0.3085, 0.2940, 0.3926, 0.3797, 0.4722, 0.4610]
CD = [0.3373, 0.3357, 0.4235, 0.4209]
CF = [0.3807, 0.3805, 0.4628, 0.4624]
CG = [0.4251, 0.4251]

BS = [0.1178, 0.0936, 0.1090, 0.1016, 0.1270, 0.1224, 0.1470, 0.1432, 0.1658,
      0.1630, 0.1823, 0.1795]
BP = [0.0874, 0.0865, 0.1099, 0.1090, 0.1319, 0.1309]
BD = [0.0960, 0.0960, 0.1193, 0.1192]
BF = [0.1078, 0.1078]

# The plot

plt.xlabel('Number of states', fontweight='bold', fontsize=15)
plt.ylabel('<v^2> / GeV', fontweight='bold', fontsize=15)

PCS = plt.scatter(np.arange(1,9), CS, s=100, color='r', marker='1')
PCP = plt.scatter(np.arange(1,7), CP, s=100, color='r', marker='2')
PCD = plt.scatter(np.arange(1,5), CD, s=100, color='r', marker='3')
PCF = plt.scatter(np.arange(1,5), CF, s=100, color='r', marker='4')
PCG = plt.scatter(np.arange(1,3), CG, s=100, color='r', marker='+')

PBS = plt.scatter(np.arange(1,13), BS, s=100, color='b', marker='1')
PBP = plt.scatter(np.arange(1,7), BP, s=100, color='b', marker='2')
PBD = plt.scatter(np.arange(1,5), BD, s=100, color='b', marker='3')
PBF = plt.scatter(np.arange(1,3), BF, s=100, color='b', marker='4')

plt.axhline(0, color='black')
plt.axis([0, 14, 0, 0.6])

leg = plt.legend((PCS, PCP, PCD, PCF, PCG),
                 ('S-States', 'P-States', 'D-States', 'F-States', 'G-States'),
                 scatterpoints=1, loc='upper right', ncol=1, fontsize=15)
plt.gca().add_artist(leg)

plt.legend((PBS, PBP, PBD, PBF), ('S-States', 'P-States', 'D-States',
                                   'F-States'), scatterpoints=1,
           loc='lower right', ncol=1, fontsize=15)

plt.show()
      
