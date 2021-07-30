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
date_str = "21-07-26_11-42-06"
u_path_str = "../data_bckp_" + date_str + "/file_SVs_circle_200_1e-16_single.dat"
u_path = os.path.join(u_path_str)
i = 0
for algorithm_str in ["Brent", "Newton"]:
    fig = plt.figure()
    fig.suptitle(r"Minima found via " + algorithm_str + " using 200 panels for operators, $c_o=20, c_i=1$")
    for num_steps in [10, 20, 50, 100]:
        roots_path_str = "../data_bckp_"+ date_str +     "/file_roots_"+ algorithm_str.lower() + "_circle_arnoldi_20_"+ str(num_steps) +     "_1e-2.dat"
        roots_path = os.path.join(roots_path_str)
        sample_size = np.loadtxt(u_path).shape[0]
        u = np.zeros((sample_size, 2))
        roots = np.zeros((sample_size, 5))
        u = np.loadtxt(u_path)
        roots = np.loadtxt(roots_path)
        roots = roots[~np.isnan(roots).any(axis=1)]
        roots_size = roots.shape[0]
        toBeRemoved = np.zeros(roots_size, dtype=bool)
        #for j in np.arange(roots_size-1):
        #    if roots[j,3]<roots[j+1,3]:
        #        toBeRemoved[j+1] = True
        #roots = roots[~toBeRemoved]
        print(algorithm_str, num_steps)
        print(roots)
        plt.subplot(2,2,(i)%4+1)
        if i%2 == 0:
            plt.ylabel("smallest singular value and roots")
        if i > 1:
            plt.xlabel("wavenumber k")
        plt.plot(u[:, 0], u[:, 1])
        plt.plot(roots[:, 1], roots[:, 3], "o", color="red")
        plt.title(r"" + str(num_steps) + " grid points on k for     root search")
        i += 1

    plt.tight_layout()
    plt.savefig("../figures/minima_search_" + algorithm_str.lower() + ".pdf")
    plt.show()
    plt.close(fig)
    #iter_path[i] = os.path.join("../data/file_iter_1e-" + str(pow(2,-(i-4))) +".dat")
#for i in range(0,3):
#    u_path[5+i] = os.path.join("../data/file_vals_1e" + str(i) +".dat")
#    roots_path[5+i] = os.path.join("../data/file_roots_1e" + str(i) +".dat")
    #iter_path[5+i] = os.path.join("../data/file_iter_1e" + str(i) +".dat")



#iter = np.zeros((sample_size,8*3))
#for i in range(8):
#    u[:,2*i:2*i+2] = np.loadtxt(u_path[i])
#    roots[:,5*i:5*i+5] = np.loadtxt(roots_path[i])
    #iter[:,3*i:3*i+3] = np.loadtxt(iter_path[i])


#for i in range(0,5):
#    if i%4==0:
#        fig.suptitle("Comparing Minima")
#    plt.subplot(2,2,(i)%4+1)
#    if i%2==0:
#    if i%4>1:
    #print(r"accuracy $=10^{-" + str(i) + "}$")
    #print(r"Minimum at: " + str(roots[:,5*i+1]) + "\n"
    #    "Average #restarts: " + str(np.mean(iter[:,3*i+1])) + "\n"
    #    "Average #matXvec: " + str(np.mean(iter[:, 3 * i + 2])))
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

