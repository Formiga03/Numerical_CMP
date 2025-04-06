import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def Reader_output(filename):
    """ Reads a .txt file with the quantity's information by inserting the its path + filename
    out putting a list of lists, in which each list is the values in a line."""
    with open(filename, "r") as file:
        lst1 = []
        for line in file:
            numbers = [float(num) for num in line.split()]
            lst1.append(numbers)
    return lst1

def IPR_sorter(lst):
    res =[lst[0], lst[1]]

    for ii in range(len(lst[1])):
        res.append([])
    
    for jj in range(2, len(lst)):
        for kk in range(len(lst[1])):
            res[kk+2].append(lst[jj][kk])
    return res

def plot_graph_grdst(lst, sets, filename, title="", style="-x"):
    W_data = [x for x in range(int(lst[0][0]))]

    for ii in sets:
        plt.title(title)
        plt_flname = "../plots/"+filename[:-3]+"_W={"
        for jj in ii:
            plt.plot(W_data, lst[jj+2], style, label="W="+str(lst[1][jj]))
            plt_flname += str(lst[1][jj]) + " "
        plt_flname = plt_flname[:-1] + "}.png" 
        plt.xlabel(r'$\left| \psi_i \right|$')
        plt.ylabel("lattice point i")
        plt.legend()       
        plt.savefig(plt_flname)
        plt.clf()

def plot_graph_ESpec(lst, filename, title="", style="-x", axis=["", ""]):
    dims = len(lst)-2
    beta_data = [2*np.pi*x*0.0001 for x in range(1600)]
    plt.title(title)
    plt_flname = "../plots/"+filename[:-3]
    for ii in tqdm(range(2,dims)):
        for jj in range(len(lst[ii])):
            plt.plot(beta_data[ii], lst[ii][jj], 'bo', markersize=0.05) 
    plt.xlabel(r'$\beta$')
    plt.ylabel("Energy Spectrum")
    plt.ylim(-1,1)
    plt.savefig(plt_flname)
    plt.clf()

def plot_graph(lst, sets, filename, title="", style="-x", axis=["", ""]):
    if lst[1][2]:
        t_data = [2**x for x in np.arange(lst[1][0],lst[1][1],lst[1][2])]
        for ii in sets:
            plt.title(title)
            plt_flname = "../plots/"+filename[:-3]+"_L={"
            for jj in ii:
                plt.semilogx(t_data, lst[jj+3], style, label="L="+str(lst[2][jj]))
                plt_flname += str(lst[2][jj]) + " "
            plt_flname = plt_flname[:-1] + "}.png" 
            plt.xlabel(axis[0])
            plt.ylabel(axis[1])
            plt.legend()          
            plt.savefig(plt_flname)
            plt.clf()

    if not lst[1][2]:
        t_data = [x for x in np.arange(lst[1][0],lst[1][1],lst[1][2])]

        for ii in sets:
            plt.title(title)
            plt_flname = "../plots/"+filename[:-3]+"_L={"
            for jj in ii:
                plt.plot(t_data, lst[jj+3], "-x", label="L="+str(lst[2][jj]))
                plt_flname += str(lst[2][jj]) + " "
            plt_flname = plt_flname[:-1] + "}.png"           
            plt.xlabel(axis[0])
            plt.ylabel(axis[1])
            plt.legend()          
            plt.savefig(plt_flname)
            plt.clf()

def plot_2graphs(lst1, lst2, sets1, sets2, filename, title="", style="-x", axis=[["",""],["", ""]]):

    if len(sets1) == len(sets2):
        
        for ii in range(len(sets1)):
            fig, axs = plt.subplots(1, 2, figsize=(10, 4)) 
            plt_flname = "../plots/"+filename[0][:-3]+"_L1={"

            for jj1 in range(len(sets1[ii])):
                if lst1[1][3]:
                    t_data1 = [2**x for x in np.arange(lst1[1][0],lst1[1][1],lst1[1][2])]
                    axs[0].semilogx(t_data1, lst1[sets1[ii][jj1]+3], style, label="L="+str(lst1[2][sets1[ii][jj1]]))
                    axs[0].set_xlabel(axis[0][0])
                    axs[0].set_ylabel(axis[0][1])
                    axs[0].legend()         
                
                if not lst1[1][3]:
                    t_data1 = [x for x in np.arange(lst1[1][0],lst1[1][1],lst1[1][2])]
                    axs[0].plot(t_data1, lst2[sets1[ii][jj1]+3], style, label="L="+str(lst1[2][sets1[ii][jj1]]))
                    axs[0].set_xlabel(axis[0][0])
                    axs[0].set_ylabel(axis[0][1])
                    axs[0].legend()         

                plt_flname += str(lst1[2][sets1[ii][jj1]]) + " "
            
            plt_flname = plt_flname[:-1] + "}_" + filename[1][:-3] + "_L2={"
            
            for jj2 in range(len(sets2[ii])):
                    
                if lst2[1][3]:
                    t_data2 = [2**x for x in np.arange(lst2[1][0],lst2[1][1],lst2[1][2])]
                    axs[1].semilogx(t_data2, lst2[sets2[ii][jj2]+3], style, label="L="+str(lst2[2][sets2[ii][jj2]]))
                    axs[1].set_xlabel(axis[1][0])
                    axs[1].set_ylabel(axis[1][1])
                    axs[1].legend()         
                
                if not lst2[1][3]:
                    t_data2 = [x for x in np.arange(lst2[1][0],lst2[1][1],lst2[1][2])]
                    axs[1].plot(t_data2, lst2[sets2[ii][jj2]+3], style, label="L="+str(lst2[2][sets2[ii][jj2]]))
                    axs[1].set_xlabel(axis[1][0])
                    axs[1].set_ylabel(axis[1][1])
                    axs[1].legend()          


                plt_flname += str(lst2[2][sets2[ii][jj2]]) + " "
            
            plt_flname = plt_flname[:-1] + "}.png"                 

            fig.suptitle(title, fontsize=10)
            fig.savefig(plt_flname)
            plt.close(fig)
   
def plot_graphs_IPR(lst, sets, filename, title="", style="-o"):
    W_data = [x for x in np.arange(lst[0][0],lst[0][1],lst[0][2])]

    for ii in sets:
        plt.title(title)
        plt_flname = "../plots/"+filename[:-3]+"_L={"
        for jj in ii:
            plt.plot(W_data, lst[jj+2], style, label="L="+str(lst[1][jj]))
            plt_flname += str(lst[1][jj]) + " "
        plt_flname = plt_flname[:-1] + "}.png" 
        plt.xlabel('perturbation amplitude W')
        plt.ylabel("IPR")          
        plt.legend()         
        plt.semilogy()
        plt.savefig(plt_flname)
        plt.clf()
        



    
