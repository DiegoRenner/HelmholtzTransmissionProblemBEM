import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
from io import StringIO

rc('text', usetex=True)
#plt.rc('font', size=10) #controls default text size
#plt.rc('axes', titlesize=10) #fontsize of the title
plt.rc('axes', labelsize=20) #fontsize of the x and y labels
plt.rc('xtick', labelsize=20) #fontsize of the x tick labels
plt.rc('ytick', labelsize=20) #fontsize of the y tick labels
plt.rc('legend', fontsize=20) #fontsize of the legend
# maybe legend size smaller?

date_str_minima = "22-04-14_22-47-54"
domains = ["circle", "square"]
domain_str = ""
u_path_str = ""
k_grids = [1,3,10]
panels_circle = [20, 50]
panels_square = [20, 48]
for k_grid in k_grids:
    for domain in domains:

        panels = []
        if domain == domains[0]:
            panels = panels_circle
        if domain == domains[1]:
            panels = panels_square

        for panel in panels:
            if domain == domains[0]:
                domain_str = "circle"
            if domain == domains[1]:
                domain_str = "square"
                panels = panels_square

            u_path = os.path.join(u_path_str)
            for acc in ["1e-4"]:
                 roots_path_str = ""
                 if domain == domains[0]:
                     roots_path_str = "../../data_bckp_" + date_str_minima + "/file_roots_seq_" + domain_str +     "_arnoldi_20_" + str(k_grid) + "_" + acc + "_" + str(panel) + ".dat"
                 if domain == domains[1]:
                     roots_path_str = "../../data_bckp_" + date_str_minima + "/file_roots_seq_" + domain_str +     "_arnoldi_20_" + str(k_grid) + "_" + acc + "_" + str(panel) + ".dat"

                 roots_path = os.path.join(roots_path_str)
                 roots_txt = open(roots_path)
                 roots_file = roots_txt.readlines()
                 function_calls = roots_file[-2][17:-1]
                 in_table = 0
                 active = []
                 num_active = []
                 no_zero = []
                 num_no_zero = []
                 zero_found = []
                 num_zero_found = []
                 grid_counter = 0
                 num_intervalls = []
                 temp_num_intervalls = 0
                 for line in roots_file:
                     if in_table==0:
                         if line.find("Grid Points")>=0:
                             in_table = 1
                     else:
                         if line == "\n":
                             in_table = 0
                             grid_counter = grid_counter + 1
                             num_intervalls.append(temp_num_intervalls)
                             temp_num_intervalls = 0
                         else:
                             intervall_type = np.loadtxt(StringIO(line))[3]
                             temp_num_intervalls = temp_num_intervalls + 1
                             if intervall_type== 0:
                                 active.append(np.loadtxt(StringIO(line))[0])
                                 num_active.append(grid_counter)
                             if intervall_type== 1:
                                 no_zero.append(np.loadtxt(StringIO(line))[0])
                                 num_no_zero.append(grid_counter)
                             if intervall_type== 2:
                                 if(np.loadtxt(StringIO(line))[0] > 0.1):
                                    zero_found.append(np.loadtxt(StringIO(line))[0])
                                    num_zero_found.append(grid_counter)

                 fig,ax = plt.subplots()

                 ax.set_xlabel("iteration")
                 ax.set_ylabel("search grid over wavenumber k")
                 ax.yaxis.get_major_locator().set_params(integer=True)
                 ax.xaxis.get_major_locator().set_params(integer=True)
                 ax.plot(num_active, active, "x", color="blue", ms=2)
                 ax.plot(num_no_zero, no_zero, "x", color="red", ms=2)
                 ax.plot(num_zero_found, zero_found, "x", color="green")
                 plt.legend(["active","no zero","zero found"], loc="lower right")
                 ax2 = ax.twinx()
                 ax2.plot(num_intervalls, color="pink")
                 plt.legend(["function calls"], loc="upper left")
                 ax2.set_ylabel("function calls")
                 plt.tight_layout()
                 plt.savefig("../../figures/minima_search_grid_" + domain_str.lower() + "_" +str(acc) + "_" +     str(panel) + "_" + function_calls + "_" + str(k_grid)  + ".pdf")
                 plt.savefig("../../figures/minima_search_grid_" + domain_str.lower() + "_" +str(acc) + "_" +  str(panel) + "_" + function_calls + "_" + str(k_grid) + ".eps")
                 plt.show()
                 plt.close(fig)