import sys
import os
sys.path.append(os.path.abspath("../src"))
from plotter import *


directory = "../outputs"  # Change this to your target directory
files = os.listdir(directory)
files.sort()

print(files) 
 
for ii in range(len(files)//2):
    rf = Reader_output("../outputs/" + files[ii]) 
    stdev = Reader_output("../outputs/" + files[len(files)//2 + ii]) 

    rf_dist = files[ii].split("_")
    stdev_dist = files[len(files)//2 + ii].split("_")
    print(rf_dist, stdev_dist)

    if rf_dist[-2] == "OffDiag":
        title= f"Return Probability (right) and Standard Deviation (left) with a perturbation following a\n {rf_dist[1]}, W={rf[0][0]} and off the diagonal of the hamitonian."

    else:
        title= f"Return Probability (right) and Standard Deviation (left) with \na perturbation following a {rf_dist[1]} and W={rf[0][0]}."

    plot_2graphs(rf, stdev, [[0,1],[2,3]], [[0,1],[2,3]], 
                [files[ii], files[len(files)//2 + ii]], 
                title=title)

"""

    plot_graph(return_funct, [[0,1], [2,3]], files[0],
                title="Return Prob. for a Gamma Distr. Perturbation with W=")
    plot_graph(data2, [[0,1],[0,2]], files[len(files)//2], 
            title="Normal Distr. Perturbation with W=")

"""