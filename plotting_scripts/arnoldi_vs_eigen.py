import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)

acc = input("Enter accurracy with which data was computed (1e[n]):")
for i in range(0,7):
    fname = os.path.join("../data/file_timings_" + str(pow(2,i)) + "_1e" + acc + ".dat")
    timings = np.loadtxt(fname)
    fname = os.path.join("../data/file_vals_arpp_" + str(pow(2,i)) + "_1e" + acc + ".dat")
    vals = np.loadtxt(fname)
    if i%8==0:
        fig = plt.figure(figsize=(640/100,960/100))
    if i==0:
        fig.suptitle(r"Runtime at each wavenumber k, accuracy $=10^{" + acc + "}$")


    X_timings = timings[:,0]
    Y1_timings = timings[:,1]
    Y2_timings = timings[:,2]
    Y3_timings = timings[:,3]
    Y_vals = vals[:,1]

    ax = plt.subplot(4,2,(i)%8+1)
    if i%2==0:
        plt.ylabel("runtime [ms]")
    if i==6:
        plt.xlabel("wavenumber k")
    plt.plot(X_timings, Y1_timings)
    plt.plot(X_timings, Y2_timings)
    if i==1:
        plt.legend(["direct", "Arnoldi"])
    ax_alt = ax.twinx()
    ax_alt.plot(X_timings,Y_vals, color="grey", alpha=0.5)
    if i%2==1:
        plt.ylabel("smallest SV")
    print("Requested " + str(pow(2,i)) + " SV.")
    print("direct average [ms]: " + str(np.mean(Y1_timings)))
    print("Arnoldi average [ms]: " + str(np.mean(Y2_timings)))
    plt.xlim(0, 10)
    plt.title(str(pow(2, i)) + " SV requested")
    plt.tight_layout()

fname = os.path.join("../data/file_timings_98_1e" + acc + ".dat")
timings = np.loadtxt(fname)
fname = os.path.join("../data/file_vals_arpp_98_1e" + acc + ".dat")
vals = np.loadtxt(fname)
X_timings = timings[:, 0]
Y1_timings = timings[:, 1]
Y2_timings = timings[:, 2]
Y3_timings = timings[:, 3]
Y_vals = vals[:, 1]

ax = plt.subplot(4, 2, 8)
plt.plot(X_timings, Y1_timings)
plt.plot(X_timings, Y2_timings)
plt.xlabel("wavenumber k")
ax_alt = ax.twinx()
ax_alt.plot(X_timings, Y_vals, color="grey", alpha=0.5)
plt.ylabel("smallest SV")
print("Requested 100 SV.")
print("direct average [ms]: " + str(np.mean(Y1_timings)))
print("Arnoldi average [ms]: " + str(np.mean(Y2_timings)))
plt.xlim(0, 10)
plt.title("100 SV requested")
plt.savefig("../figures/arnoldi_time_1e" + acc + ".pdf")
plt.show()
plt.close(fig)

for i in range(0,7):
    fname = os.path.join("../data/file_iter_" + str(pow(2,i)) + "_1e" + acc + ".dat")
    iter  = np.loadtxt(fname)
    fname = os.path.join("../data/file_vals_arpp_" + str(pow(2,i)) + "_1e" + acc + ".dat")
    vals = np.loadtxt(fname)
    if i==0:
        fig = plt.figure(figsize=(640 / 100, 960 / 100))
        fig.suptitle(r"Iterations at each wavenumber k, accuracy $=10^{" + acc + "}$")

    X_iter = iter[:,0]
    Y1_iter = iter[:,1]
    Y2_iter = iter[:,2]
    Y_vals = vals[:,1]

    ax = plt.subplot(4,2,i+1)
    if i%2==0:
        plt.ylabel("iterations")
    if i==6:
        plt.xlabel("wavenumber k")
    plt.plot(X_iter, Y1_iter)
    plt.plot(X_iter, Y2_iter)
    if i==1:
        plt.legend(["restarts", "matrixXvector"])
    ax_alt = ax.twinx()
    ax_alt.plot(X_timings,Y_vals,color="grey", alpha=0.5)
    if i%2==1:
        plt.ylabel("smallest SV")
    print("Requested " + str(pow(2,i)) + " SV.")
    print("restarts average iterations: " + str(np.mean(Y1_iter)))
    print("matXvec average iteration: " + str(np.mean(Y2_iter)))
    plt.xlim(0, 10)
    plt.title(str(pow(2, i)) + " SV requested")
    plt.tight_layout()

fname = os.path.join("../data/file_iter_98_1e" + acc + ".dat")
iter = np.loadtxt(fname)
fname = os.path.join("../data/file_vals_arpp_98_1e" + acc + ".dat")
vals = np.loadtxt(fname)
X_iter = iter[:, 0]
Y1_iter = iter[:, 1]
Y2_iter = iter[:, 2]
Y_vals = vals[:, 1]

ax = plt.subplot(4, 2, 8)
plt.plot(X_iter, Y1_iter)
plt.plot(X_iter, Y2_iter)
plt.xlabel("wavenumber k")
ax_alt = ax.twinx()
ax_alt.plot(X_timings, Y_vals, color="grey", alpha=0.5)
plt.ylabel("smallest SV")
print("Requested 100 SV.")
print("restarts average: " + str(np.mean(Y1_iter)))
print("matXvec average: " + str(np.mean(Y2_iter)))
plt.xlim(0, 10)
plt.title("100 SV requested")
plt.savefig("../figures/arnoldi_iter_1e" + acc + ".pdf")
plt.show()
plt.close(fig)

