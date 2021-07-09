import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)

#u_path = [None]*8
#roots_path = [None]*8
#iter_path = [None]*8
#for i in range(0,5):
    #u_path[i] = os.path.join("../data_bckp_21-06-09-20-19/file_vals_1e-" + str(pow(2,-(i-4))) +".dat")
u_path = os.path.join("../data_bckp_21-07-09_11-15-44/file_SVs_circle_50_1e-16_single.dat")
roots_path = os.path.join("../data_bckp_21-06-09-20-19/file_roots_brent_circle_arnoldi_20_100_1e-8.dat")
    #iter_path[i] = os.path.join("../data/file_iter_1e-" + str(pow(2,-(i-4))) +".dat")
#for i in range(0,3):
#    u_path[5+i] = os.path.join("../data/file_vals_1e" + str(i) +".dat")
#    roots_path[5+i] = os.path.join("../data/file_roots_1e" + str(i) +".dat")
    #iter_path[5+i] = os.path.join("../data/file_iter_1e" + str(i) +".dat")



sample_size = np.loadtxt(u_path).shape[0]
u = np.zeros((sample_size,2))
print(u[:,0])
roots = np.zeros((sample_size,5))
#iter = np.zeros((sample_size,8*3))
#for i in range(8):
#    u[:,2*i:2*i+2] = np.loadtxt(u_path[i])
u = np.loadtxt(u_path)
#    roots[:,5*i:5*i+5] = np.loadtxt(roots_path[i])
roots = np.loadtxt(roots_path)
    #iter[:,3*i:3*i+3] = np.loadtxt(iter_path[i])

roots=roots[~np.isnan(roots).any(axis=1)]

#for i in range(0,5):
#    if i%4==0:
#        fig = plt.figure()
#        fig.suptitle("Comparing Minima")
#    plt.subplot(2,2,(i)%4+1)
#    if i%2==0:
#        plt.ylabel("smallest singular value")
#    if i%4>1:
#        plt.xlabel("wavenumber k")
plt.plot(u[:,0], u[:,1])
print(u[:,0])
plt.plot(roots[:,1], roots[:,3], "o",color="red")

    plt.title(r"Roots found ")
    #print(r"accuracy $=10^{-" + str(i) + "}$")
    #print(r"Minimum at: " + str(roots[:,5*i+1]) + "\n"
    #    "Average #restarts: " + str(np.mean(iter[:,3*i+1])) + "\n"
    #    "Average #matXvec: " + str(np.mean(iter[:, 3 * i + 2])))
plt.tight_layout()
#    if i==3:
#        plt.show()
#        plt.close(fig)
#
#for i in range(0,3):
#    if (i+5)%4==0:
#        fig = plt.figure()
#        fig.suptitle("Comparing Minima")
#    plt.subplot(2,2,(i+5)%4+1)
#    if (i+5)%2==0:
#        plt.ylabel("smallest singular value")
#    if (i+5)%4>1:
#        plt.xlabel("wavenumber k")
#    plt.plot(u[:,2*(i+5)], u[:,2*(i+5)+1])
#    plt.plot(roots[:,5*(i+5)+1], roots[:,5*(i+5)+3], "o",color="red")
#    plt.title(r"accuracy $=10^{" + str(i) + "}$")
#    print(r"accuracy $=10^{" + str(i) + "}$")
#    print("Minimum at: " + str(roots[:,5*(i+5)+1]) + "\n"
#        "Average #restarts: " + str(np.mean(iter[:,3*(i+5)+1])) + "\n"
#        "Average #matXvec: " + str(np.mean(iter[:, 3 * (i+5) + 2])))
#    plt.tight_layout()
#    if (i+5)%4==3:
#        plt.show()
#        plt.close(fig)
#
plt.show()
#plt.close(fig)

