import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)


l = 1e-16 #input("Enter accurracy for arnoldi algorithm:")
X = np.zeros(6)
k1_values = np.zeros((6,8))
k2_values = np.zeros((6,8))
k3_values = np.zeros((6,8))
k1_fourier = np.array([0.60669079, 0.96913765, 0.99047821, 0.99570781, 1.00430351, 1.00954957,
 1.03264829, 1.73330956], np.float64)
k2_fourier = np.array([0.49797674, 0.67860537, 0.75321262, 0.78150035, 0.81515499, 0.85032776,
 0.87243619, 0.88814115],np.float64)
k3_fourier = np.array([0.13278006, 0.41964322, 0.43699315, 0.44207863, 0.48421406, 0.49793782,
 0.51548511, 0.52753431], np.float64)
values = np.zeros((3,8))
i = 0
for k in [50,100,200,400,800,1600]:
    #fig = plt.figure()
    #fig.suptitle("Disk-shaped scatterer: inner index = 3.0, outer index = 1.0\n "
    #        r"\#panels used for operator approximation: " + str(k) +
    #        r", accurracy for arnoldi algorithm: " + str(l))
    fname = os.path.join("../data/file_SVs_" + str(k)  + "_" + str(l) + ".dat")
    values = np.loadtxt(fname)
    values = values[values[:,0].argsort()]

    X[i] = 2*np.pi/k
    k1_values[i,:] = abs(values[0,1:9] - k1_fourier)
    k2_values[i,:] = abs(values[1,1:9] - k2_fourier)
    k3_values[i,:] = abs(values[2,1:9] - k3_fourier)
    print(k1_values)
    i+=1
    #X_residuals = residuals[:, 0]
    #Y1_residuals = residuals[:, 1:7]
    #
    #plt.ylabel("singular values")
    #plt.xlabel("wave number k")
    #plt.plot(X_residuals, Y1_residuals)
    #plt.tight_layout()
    #plt.savefig("../figures/SVs_" + str(k) + "_"+ str(l) + ".pdf")
    #plt.show()
    #plt.close(fig)

fig,ax = plt.subplots()
fig.suptitle("Disk-shaped scatterer: inner index = 3.0, outer index = 1.0\n "
             r"k = " + str(values[0,0]) +
             r", arnoldi algorithm accurracy: " + str(l))
plt.loglog(X,k1_values[:,0])
ax.invert_xaxis()
plt.xlabel("meshwidth h [log]")
plt.ylabel(r"$|\lambda_{0,fourier}-\lambda_{0,BEM}|$ [log]")
plt.tight_layout()
plt.savefig("../figures/SVs_" + str(values[0,0]) + "_" + str(l) + ".pdf")
plt.show()
plt.close(fig)

fig,ax = plt.subplots()
fig.suptitle("Disk-shaped scatterer: inner index = 3.0, outer index = 1.0\n "
             r"k = " + str(values[1,0]) +
             r", arnoldi algorithm accurracy: " + str(l))
plt.loglog(X,k2_values[:,0])
ax.invert_xaxis()
plt.xlabel("meshwidth h [log]")
plt.ylabel(r"$|\lambda_{0,fourier}-\lambda_{0,BEM}|$ [log]")
plt.tight_layout()
plt.savefig("../figures/SVs_" + str(values[1,0]) + "_" + str(l) + ".pdf")
plt.show()
plt.close(fig)

fig,ax = plt.subplots()
fig.suptitle("Disk-shaped scatterer: inner index = 3.0, outer index = 1.0\n "
             r"k = " + str(values[2,0]) +
             r", arnoldi algorithm accurracy: " + str(l))
print(k3_values[:,0])
plt.loglog(X,k3_values[:,0])
ax.invert_xaxis()
plt.xlabel("meshwidth h [log]")
plt.ylabel(r"$|\lambda_{0,fourier}-\lambda_{0,BEM}|$ [log]")
plt.tight_layout()
plt.savefig("../figures/SVs_" + str(values[2,0]) + "_" + str(l) + ".pdf")
plt.show()
plt.close(fig)
