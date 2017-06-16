import matplotlib.pyplot as plt
import numpy as np

# Plots the potential function and shows where the quarks lie with respect to it.

# Read from command line q= 1 or 2, 1 for Charmonium and 2 for Bottomonium
q = int(raw_input("1 for Charmonium or 2 for Bottomonium:"))
k = int(raw_input("Give value for k:"))
V = np.zeros(k-1)
x = np.arange(1,k)

hc = 197.327

# Create the potential function and stuff
if q==1:

    for i in range(1,k):
        V[i-1] = -4*0.54174*hc/(3*i)+0.15026*i/hc+0.32367

    S = [395, 449, 878, 1261, 1590, 910, 1284, 1607]
    P = [707, 720, 1114, 1459]
    D = [936, 938, 1127, 1471, 1308, 1311]
    F = [1130, 1130, 1478, 1478]

    VS = np.zeros(8)
    VP = np.zeros(4)
    VD = np.zeros(6)
    VF = np.zeros(4)
    
elif q==2:

    for i in range(1,k):
        V[i-1] = -4*0.33296*hc/(3*i)+0.23582*i/hc+0.98732

    S = [226, 245, 500, 723, 916, 1089, 1232, 515, 735, 925, 1098, 1238]
    P = [406, 407, 642, 843, 644, 845]
    D = [537, 537, 753, 753]
    F = [650, 650]

    VS = np.zeros(12)
    VP = np.zeros(6)
    VD = np.zeros(4)
    VF = np.zeros(2)

# Plot

plt.xlabel('Radius / am', fontweight='bold', fontsize=15)
plt.ylabel('V(r) / GeV', fontweight='bold', fontsize=15)
plt.plot(x, V, 'b-', linewidth=2.0)
if q==1:

    for i in range(0,8):
        VS[i] = V[S[i]-1]

    for i in range(0,4):
        VP[i] = V[P[i]-1]
        VF[i] = V[F[i]-1]

    for i in range(0,6):
        VD[i] = V[D[i]-1]

elif q==2:

    for i in range(0,12):
        VS[i] = V[S[i]-1]

    for i in range(0,6):
        VP[i] = V[P[i]-1]

    for i in range(0,4):
        VD[i] = V[D[i]-1]

    for i in range(0,2):
        VF[i] = V[F[i]-1]
        
PS = plt.scatter(S, VS, s=70, color='r', marker='D')
PP = plt.scatter(P, VP, s=70, color='g', marker='o')
PD = plt.scatter(D, VD, s=70, color='m', marker='s')
PF = plt.scatter(F, VF, s=70, color='k', marker='*')
plt.axhline(0, color='black')
plt.axis([0, k, -5, 5])
plt.legend((PS, PP, PD, PF), ('S-States', 'P-States', 'D-States', 'F-States'),
           scatterpoints=1, loc='lower right', ncol=1, fontsize=15)

plt.show()
    
        
    
