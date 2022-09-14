import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
from io import StringIO
rc('text', usetex=True)
rc('font', size=20)

date_str_minima = "22-03-31_21-05-33"
domains = ["circle", "square"]
domain_str = ""
u_path_str = ""
k_grids = [1,3,10]
panels_circle = []
panels_square = []
for k_grid in k_grids:
    if k_grid == k_grids[0]:
        panels_circle = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
        panels_square = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192]
    if k_grid == k_grids[1]:
        panels_circle = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200]
        panels_square = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 136, 144, 152, 160, 168, 176, 184, 192]
    if k_grid == k_grids[2]:
        panels_circle = [10, 20, 30, 40, 50, 60, 70, 80]
        panels_square = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112]
    for domain in domains:
        total_function_calls = 0
        panels = []
        if domain == domains[0]:
            panels = panels_circle
        if domain == domains[1]:
            panels = panels_square

        fig = plt.figure()
        plt.ylabel("roots found")
        plt.xlabel("\#panels")

        for i, panel in enumerate(panels):

            if domain == domains[0]:
                domain_str = "circle"
            if domain == domains[1]:
                domain_str = "square"
                panels = panels_square

            for acc in ["1e-4"]:
                roots_path_str = ""
                if domain == domains[0]:
                    roots_path_str = "../../data_bckp_" + date_str_minima +     "/file_roots_seq_" + domain_str +     "_arnoldi_20_" + str(k_grid) + "_" + acc + "_" + str(panel) + ".dat"
                if domain == domains[1]:
                    roots_path_str = "../../data_bckp_" + date_str_minima +     "/file_roots_seq_" + domain_str +     "_arnoldi_20_" + str(k_grid) + "_" + acc + "_" + str(panel) + ".dat"


                roots_path = os.path.join(roots_path_str)
                roots_txt = open(roots_path)
                roots_file = roots_txt.readlines()
                function_calls = roots_file[-2][17:-1]
                total_function_calls = total_function_calls + int(function_calls)
                last_line = StringIO(roots_file[-1][13:])
                roots = np.loadtxt(last_line)
                roots_size = roots.shape[0]
                values = roots_size*[panel]
                plt.plot(values[1:], roots[1:], "x", color="red", )
                plt.tight_layout()

        average_function_calls = round(total_function_calls / len(panels))
        plt.tight_layout()
        plt.savefig("../../figures/minima_search_conv_" + domain_str.lower() + "_" +str(acc) + "_" + str(average_function_calls) + "_" + str(k_grid)  + ".pdf")
        plt.savefig("../../figures/minima_search_conv_" + domain_str.lower() + "_" +str(acc) + "_" + str(average_function_calls) + "_" + str(k_grid) + ".eps")
        plt.show()
        plt.close(fig)