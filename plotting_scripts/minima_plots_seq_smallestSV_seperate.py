import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
from io import StringIO
rc('text', usetex=True)

date_str_minima = "22-02-25_12-01-54"
date_str_sv = "21-07-26_11-42-06"
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
        if domain == domains[0]:
            #fig.suptitle(r"Minima found using " + str(panel) + " panels for operators" + "\n"
            #            r"circle domain, $r=0.25, c_i=20, c_o=1$")
            domain_str = "circle"
            u_path_str = "../data_bckp_" + date_str_sv + "/file_SVs_" + domain_str + "_200_1e-16_single.dat"
        if domain == domains[1]:
            #fig.suptitle(r"Minima found using " + str(panel) + " panels for operators" + "\n"
            #            r"square domain, $s=0.25, c_i=20, c_o=1$")
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

             roots_path = os.path.join(roots_path_str)
             roots_txt = open(roots_path)
             roots_file = roots_txt.readlines()
             function_calls = roots_file[-2][17:-1]
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
            #plt.subplot(2,2,(i)%4+1)
             fig = plt.figure()
             plt.ylabel("smallest singular value and roots")
             plt.xlabel("wavenumber k")
             plt.plot(u[1:, 0], u[1:, 1])
             plt.plot(roots[1:], values[1:], "x", color="red")
             plt.savefig("../figures/minima_search_" + domain_str.lower() + "_" +str(acc) + "_" + str(panel) + "_" + function_calls + ".pdf")
             plt.savefig("../figures/minima_search_" + domain_str.lower() + "_" +str(acc) + "_" + str(panel) + "_" + function_calls + ".eps")
             plt.show()
             plt.close(fig)