import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
from io import StringIO
rc('text', usetex=True)

#u_path = [None]*8
#roots_path = [None]*8
#iter_path = [None]*8
#for i in range(0,5):
    #u_path[i] = os.path.join("../data_bckp_21-06-09-20-19/file_vals_1e-" + str(pow(2,-(i-4))) +".dat")
date_str_minima = "22-02-25_12-01-54"
date_str_sv = "21-07-26_11-42-06"
i = 0
domains = ["circle", "square"]
domain_str = ""
u_path_str = ""
panels_circle = [20, 50]
panels_square = [20, 48]
for domain in domains:

    panels = []
    if domain == domains[0]:
        panels = panels_circle
    if domain == domains[1]:
        panels = panels_square

    for panel in panels:
        fig = plt.figure()
        if domain == domains[0]:
            fig.suptitle(r"Minima found using " + str(panel) + " panels for operators" + "\n"
                        r"circle domain, $r=0.25, c_i=20, c_o=1$")
            domain_str = "circle"
            u_path_str = "../data_bckp_" + date_str_sv + "/file_SVs_" + domain_str + "_200_1e-16_single.dat"
        if domain == domains[1]:
            print("test")
            fig.suptitle(r"Minima found using " + str(panel) + " panels for operators" + "\n"
                        r"square domain, $s=0.25, c_i=20, c_o=1$")
            domain_str = "square"
            u_path_str = "../data_bckp_" + date_str_sv + "/file_SVs_" + domain_str + "_192_1e-16_single.dat"
            panels = panels_square

        u_path = os.path.join(u_path_str)
        for acc in ["1e-2", "1e-4", "1e-8", "1e-12"]:
             roots_path_str = ""
             if domain == domains[0]:
                 roots_path_str = "../data_bckp_" + date_str_minima + "/file_roots_seq_" + domain_str + "_arnoldi_20_10_" + acc + "_" + str(panel) + ".dat"
             if domain == domains[1]:
                 roots_path_str = "../data_bckp_" + date_str_minima + "/file_roots_seq_" + domain_str + "_arnoldi_20_10_" + acc + "_" + str(panel) + ".dat"

             print(roots_path_str)
             roots_path = os.path.join(roots_path_str)
             print(roots_path)
             print(u_path)
             roots_txt = open(roots_path)
             roots_file = roots_txt.readlines()
             function_calls = roots_file[-2][17:]
             last_line = StringIO(roots_file[-1][13:])
             u = np.loadtxt(u_path)
             roots = np.loadtxt(last_line) #np.array(list(last_line), dtype="float32")
             roots_size = roots.shape[0]
             values = np.zeros((roots_size))
             j = 0
             for root in roots:
                 idx = np.abs(u[:,0]-root).argmin()
                 values[j] = u[idx, 1]
                 j += 1
             #for j in np.arange(roots_size-1):
             #    if roots[j,3]<roots[j+1,3]:
             #        toBeRemoved[j+1] = True
             #roots = roots[~toBeRemoved]
             print(domain_str, acc)
             print(roots)
             plt.subplot(2,2,(i)%4+1)
             if i%2 == 0:
                 plt.ylabel("smallest singular value and roots")
             if i > 1:
                 plt.xlabel("wavenumber k")
             plt.plot(u[1:, 0], u[1:, 1])
             plt.plot(roots[1:], values[1:], "x", color="red")
             #plt.title(r"" + str(num_steps) + " grid points on k for     root search")
             plt.title(r"accuracy arnoldi alg: " + acc + ", fct calls: " + function_calls)
             i += 1

        plt.tight_layout()
        plt.savefig("../figures/minima_search_" + domain_str.lower() + "_" + str(panel) + ".pdf")
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

