import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc
rc('text', usetex=True)
#rc('font', size=15)


#date_strs=["_bckp_21-11-04_23-04-33",
#           "_bckp_21-11-16_20-11-44",
#           "_bckp_21-11-23_21-06-45"]
date_strs=[ "_bckp_21-11-16_20-11-44"]
for date_str in date_strs:
    #accs = ["1e0", "1e-16"];  # input("Enter accurracy for arnoldi algorithm:")
    accs = ["1e-16"];  # input("Enter accurracy for arnoldi algorithm:")
    example_str = "ex2"
    for acc in accs:
        #domain_strs = ["circle", "square"]
        domain_strs = ["square"]
        for domain_str in domain_strs:
            fig = plt.figure()
            #if domain_str == domain_strs[0]:
            #    panels = [50,100,200,400]
                #if date_str == date_strs[0]:
                #    example_str = "ex1"
                #    fig.suptitle( r"disk-shaped scatterer: $r = 0.25, c_i = 5.0, #c_o = 1.0$ " + "\n" +
                #                  r"accurracy for arnoldi algorithm: " + str(acc))
                #if date_str == date_strs[1]:
                #    example_str = "ex2"
                #    fig.suptitle( r"disk-shaped scatterer: $r = 1.0, c_i = 3.0, c_o# = 1.0$ " + "\n" +
                #                  r"accurracy for arnoldi algorithm: " + str(acc))
                #if date_str == date_strs[2]:
                #    example_str = "ex3"
                #    fig.suptitle( r"disk-shaped scatterer: $r = 0.25, c_i = 25.0, #c_o = 1.0$ " + "\n" +
                #                  r"accurracy for arnoldi algorithm: " + str(acc))
            #if domain_str == domain_strs[1]:
            panels = [48, 96, 192, 384]
                #if date_str == date_strs[0]:
                #    example_str = "ex1"
                #    fig.suptitle( r"square-shaped scatterer: $s = 0.25, c_i = 5.0, #c_o = 1.0$ " + "\n" +
                #                  r"accurracy for arnoldi algorithm: " + str(acc))
                #if date_str == date_strs[1]:
                #    example_str = "ex2"
                #    fig.suptitle( r"square-shaped scatterer: $s = 1.0, c_i = 3.0, #c_o = 1.0$ " + "\n" +
                #                  r"accurracy for arnoldi algorithm: " + str(acc))
                #if date_str == date_strs[2]:
                #    example_str = "ex3"
                #    fig.suptitle( r"square-shaped scatterer: $s = 0.25, c_i = 25.0, c_o = 1.0$ " + "\n" +
                #                  r"accurracy for arnkoldi algorithm: " + str(acc))

            for k in panels:
                fname = os.path.join("../../data" + date_str + "/file_SVs_" + domain_str + "_" + str(k) + "_" + str(acc) + "_single.dat")
                residuals = np.loadtxt(fname)

                X_residuals = residuals[:, 0]
                Y1_residuals = residuals[:, 1:9]

                plot_index = int(np.log(k/48)/np.log(2))+1
                plt.subplot(2, 2, plot_index)
                plt.title(r"\#panels = " + str(k))
                if (plot_index-1)%2 == 0:
                    plt.ylabel("singular values")
                if plot_index > 2:
                    plt.xlabel("wave number k")
                ax = plt.plot(X_residuals, np.flip(Y1_residuals,axis=1), linewidth=0.5, alpha=1)
                ax[0].set_color("grey")
                ax[1].set_color("pink")
                ax[2].set_color("brown")
                ax[3].set_color("purple")
                ax[4].set_color("red")
                ax[5].set_color("green")
                ax[6].set_color("orange")
                ax[7].set_color("blue")
            plt.tight_layout()
            plt.savefig("../../figures/SVs_" + domain_str + "_" + str(acc) + "_" + example_str + "_new.pdf")
            plt.show()
            plt.close(fig)
