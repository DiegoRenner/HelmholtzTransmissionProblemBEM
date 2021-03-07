import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)


l = input("Enter accurracy for arnoldi algorithm:")
for k in [50,100,200]:
    fig = plt.figure()
    fig.suptitle("Disk-shaped scatterer: inner index = 3.0, outer index = 1.0\n "
            r"\#panels used for operator approximation: " + str(k) +
            r", accurracy for arnoldi algorithm: " + str(l))
    fname = os.path.join("../data/file_SVs_" + str(k)  + "_" + str(l) + ".dat")
    residuals = np.loadtxt(fname)
    
    X_residuals = residuals[:, 0]
    Y1_residuals = residuals[:, 1:7]
    
    plt.ylabel("singular values")
    plt.xlabel("wave number k")
    plt.plot(X_residuals, Y1_residuals)
    plt.tight_layout()
    plt.savefig("../figures/SVs_" + str(k) + "_"+ str(l) + ".pdf")
    plt.show()
    plt.close(fig)
    

