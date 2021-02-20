import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)

l = input("Enter number of coefficients in series expansion for solution:")
#fig = plt.figure(figsize=(640 / 100, 960 / 100))
fig = plt.figure()
fig.suptitle("Convergence analysis "
             "- computed vs. exact solution on circular domain,\n "
             r"Incoming wave $ u_{inc}(r,\phi) =\sum_{m=-\l}^{\l} "
             r"a_m J_m(k\sqrt{n_o}r)e^{im\phi} $, l = " + str(l))
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
    plt.title("order of GaussQR: " + str(8+i))
    slope, intercept = np.polyfit(np.log(X_residuals), np.log(Y1_residuals), 1)
    #plt.text(0,0,str(slope))
    plt.tight_layout()
    #print(ax.xlim())
    plt.text(0.65, 0.75,"convergence rate: "
             + str(round(slope,4)
                   if abs(round(slope,4))>np.finfo(float).eps
                   else 0.0),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
plt.savefig("../figures/order_" + str(l) + "_comp.pdf")
plt.show()
plt.close(fig)

fig = plt.figure()

fig.suptitle("Convergence analysis "
             "- closest approximation in p.w. cont. FEM-space vs. exact solution,\n "
             r"Incoming wave $ u_{inc}(r,\phi) =\sum_{m=-\l}^{\l} "
             r"a_m J_m(k\sqrt{n_o}r)e^{im\phi} $, l = " + str(l))
for i in range(0,6):
    fname = os.path.join("../data/file_order_" + str(8+i) + "_" + str(l) + ".dat")
    residuals = np.loadtxt(fname)

    X_residuals = residuals[:, 0]
    Y1_residuals = residuals[:, 2]

    ax = plt.subplot(3,2,(i)+1)
    if i%2==0:
        plt.ylabel("error [log]")
    if i>3:
        plt.xlabel("meshwidth h [log]")
    plt.loglog(X_residuals, Y1_residuals)
    ax.invert_xaxis()
    plt.title("order of GaussQR: " + str(8+i))
    slope, intercept = np.polyfit(np.log(X_residuals), np.log(Y1_residuals), 1)
    #plt.text(0,0,str(slope))
    plt.tight_layout()
    #print(ax.xlim())
    plt.text(0.65, 0.75,"convergence rate: "
             + str(round(slope,4)
                   if abs(round(slope,4))>np.finfo(float).eps
                   else 0.0),
             horizontalalignment='center',
             verticalalignment='center',
             transform = ax.transAxes)
plt.savefig("../figures/order_" + str(l) + "_closest.pdf")
plt.show()
plt.close(fig)
