import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)

l = input("Enter number of coefficients in series expansion for solution:")
#fig = plt.figure(figsize=(640 / 100, 960 / 100))
fig = plt.figure()
fig.suptitle("Convergence analysis on circular domain, \#coefficients = " + str(l))
for i in range(0,6):
    fname = os.path.join("../data/file_order_" + str(8+i) + "_" + str(l) + ".dat")
    residuals = np.loadtxt(fname)


    X_residuals = residuals[:, 0]
    Y1_residuals = residuals[:, 1]

    ax = plt.subplot(3,2,(i)+1)
    if i%2==0:
        plt.ylabel("error [log]")
    if i>3:
        plt.xlabel("meshwidth h [log]")
    plt.loglog(X_residuals, Y1_residuals)
    ax.invert_xaxis()
    plt.title("order = " + str(8+i))
    slope, intercept = np.polyfit(np.log(X_residuals), np.log(Y1_residuals), 1)
    #plt.text(0,0,str(slope))
    plt.tight_layout()
    #print(ax.xlim())
    plt.text(0.7, 0.7,"convergence rate: \n" + str(round(slope,4)),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
plt.savefig("../figures/order_ " + str(l) + ".pdf")
plt.show()
plt.close(fig)

