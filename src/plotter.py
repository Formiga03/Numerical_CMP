import numpy as np
import matplotlib.pyplot as plt

def Reader_output(filename):
    with open(filename, "r") as file:
        lst1 = []
        for line in file:
            numbers = [float(num) for num in line.split()]
            lst1.append(numbers)
    return lst1

def plot_graph(lst, sets, filename, title="", style="-x"):

    if lst[1][3]:
        t_data = [2**x for x in np.arange(lst[1][0],lst[1][1],lst[1][2])]
        for ii in sets:
            plt.title(title)
            plt_flname = "../plots/"+filename[:-3]+"_L={"
            for jj in ii:
                plt.semilogx(t_data, lst[jj+3], style, label="L="+str(lst[2][jj]))
                plt_flname += str(lst[2][jj]) + " "
            plt_flname = plt_flname[:-1] + "}.png" 
            plt.legend()          
            plt.savefig(plt_flname)
            plt.clf()

    if not lst[1][3]:
        t_data = [x for x in np.arange(lst[1][0],lst[1][1],lst[1][2])]

        for ii in sets:
            plt.title(title)
            plt_flname = "../plots/"+filename[:-3]+"_L={"
            for jj in ii:
                plt.plot(t_data, lst[jj+3], "-x", label="L="+str(lst[2][jj]))
                plt_flname += str(lst[2][jj]) + " "
            plt_flname = plt_flname[:-1] + "}.png"           
            plt.legend()          
            plt.savefig(plt_flname)

def plot_2graphs(lst1, lst2, sets1, sets2, filename, title="", style="-x"):

    if len(sets1) == len(sets2):
        
        for ii in range(len(sets1)):
            print("ii=" +str(ii))
            print(sets1[ii], sets2[ii])
            fig, axs = plt.subplots(1, 2, figsize=(10, 4)) 
            plt_flname = "../plots/"+filename[0][:-3]+"_L1={"

            for jj1 in range(len(sets1[ii])):
                if lst1[1][3]:
                    t_data1 = [2**x for x in np.arange(lst1[1][0],lst1[1][1],lst1[1][2])]
                    axs[0].semilogx(t_data1, lst1[sets1[ii][jj1]+3], style, label="L="+str(lst1[2][sets1[ii][jj1]]))
                    axs[0].legend()         
                
                if not lst1[1][3]:
                    t_data1 = [x for x in np.arange(lst1[1][0],lst1[1][1],lst1[1][2])]
                    axs[0].plot(t_data1, lst2[sets1[ii][jj1]+3], style, label="L="+str(lst1[2][sets1[ii][jj1]]))
                    axs[0].legend()         

                plt_flname += str(lst1[2][sets1[ii][jj1]]) + " "
            
            plt_flname = plt_flname[:-1] + "}_" + filename[1][:-3] + "_L2={"
            
            for jj2 in range(len(sets2[ii])):
                    
                if lst2[1][3]:
                    t_data2 = [2**x for x in np.arange(lst2[1][0],lst2[1][1],lst2[1][2])]
                    axs[1].semilogx(t_data2, lst2[sets2[ii][jj2]+3], style, label="L="+str(lst2[2][sets2[ii][jj2]]))
                    axs[1].legend()         
                
                if not lst2[1][3]:
                    t_data2 = [x for x in np.arange(lst2[1][0],lst2[1][1],lst2[1][2])]
                    axs[1].plot(t_data2, lst2[sets2[ii][jj2]+3], style, label="L="+str(lst2[2][sets2[ii][jj2]]))
                    axs[1].legend()          


                plt_flname += str(lst2[2][sets2[ii][jj2]]) + " "
            
            plt_flname = plt_flname[:-1] + "}.png"                 

            fig.suptitle(title, fontsize=10)
            fig.savefig(plt_flname)
            plt.close(fig)

            



    
