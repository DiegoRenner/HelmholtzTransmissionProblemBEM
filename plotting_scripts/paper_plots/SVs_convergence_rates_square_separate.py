import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import rc

rc('text', usetex=True)
rc('font', size=20)
#rc('font', size=15)


date_strs=[ "_bckp_22-05-06_14-03-36"]
#date_strs=[ "_bckp_21-06-09-20-19"]
for date_str in date_strs:
    accs = ["1e-16"];  # input("Enter accurracy for arnoldi algorithm:")
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
            panels = [48, 96, 192, 384, 768, 1536]

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

            first_SV = 1
            n_SVs = 4
            last_SV = first_SV + n_SVs
            data = np.zeros((4,n_SVs,len(panels)))
            wave_number = [0.5, 1.0, 1.5, 2.0]

            for i,n in enumerate(panels):
                fname = os.path.join("../../data" + date_str     + "/file_SVs_" + domain_str + "_" +     str(n) + "_" + str(acc) +     "_single.dat")
                #fname = os.path.join("../../data" + date_str     + "/file_SVs_" + domain_str + "_" +     str(n) + "_" + str(acc) +     ".dat")
                data_raw = np.loadtxt(fname)
                print(fname)
                print(data_raw)
                data[:,:,i] = data_raw[:,first_SV:last_SV]
                #data[:, :, 0:5] = np.abs(data[:, :, 0:5] - data[:, :, 5])

            for i,k in enumerate(wave_number):
                plt.title(r"wave number $k=" + str(k) + "$")
                plt.xlabel(r"\#panels")
                plt.ylim([0.00001, 0.1])
                plt.xlim([40, 1000])
                for j in range(0,4):
                    data[i, j, 0:5] = np.abs(data[i, j, 0:5] - data[i, j, 5])
                    plot_index = i
                    #plt.subplot(2, 2, plot_index+1)
                    #if (plot_index)%2 == 0:
                    plt.ylabel("Error in singular value")
                    ax = plt.loglog(panels[0:5], data[i,j,0:5], "+", linestyle="-", label=r"$\sigma_" + str(j+1) + "$")
                    #ax = plt.plot(X_residuals, np.flip(        Y1_residuals,axis=1), linewidth=0.5,         alpha=1)
                    #ax[0].set_color("grey")
                    #ax[1].set_color("pink")
                    #ax[2].set_color("brown")
                    #ax[3].set_color("purple")
                    #ax[4].set_color("red")
                    #ax[5].set_color("green")
                    #ax[6].set_color("orange")
                    #ax[7].set_color("blue")
                plt.legend()
                plt.tight_layout()
                plt.savefig("../../figures/"+ domain_str + "_SVsCvg_" + str(k) + ".pdf")
                plt.show()
                plt.close(fig)
