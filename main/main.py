import sys
import os

sys.path.append(os.path.abspath("../src"))
from plotter import *

from tqdm import tqdm

directory = "../outputs" 
files = os.listdir(directory)
files.sort()

rf_sd = []
imbalance = []
peri = []
perialp = []
S = []
IPR = []
IPR_W = []
grdst = []
grdst_PBC = []
ESpec = []

for xx in files:
    if "Imbalance" in xx:
        imbalance.append(xx)

    if  "StandDev" in xx:
        rf_sd.append(xx)

    if "ReturnProb" in xx:
        rf_sd.append(xx)

    if "QuasiPeriodicPert" in xx:
        if "IPR" in xx:
            pass
        if "GroundState" in xx:
            if "PBC" in xx:
                grdst_PBC.append(xx)
            if not "PBC" in xx:
                grdst.append(xx)
        if "EnergySpectrum" in xx:
            ESpec.append(xx)
        elif "QuasiPeriodicPert_withAlpha_" in xx:
            perialp.append(xx)
        else:
            peri.append(xx)

    if "EntenglEntrop" in xx:
        S.append(xx)

    if "IPR" in xx:
        if "Wvs" in xx:
            IPR_W.append(xx)
        else:
            IPR.append(xx)


for ii in tqdm(ESpec):
    ESpec_dt = Reader_output("../outputs/"+ii)
    plot_graph_ESpec(ESpec_dt, ii, style="o")

if rf_sd:
    print("Return Probability and Standard Deviation for a lattice with a radom perturbation field:")
    for ii in tqdm(range(len(rf_sd)//2)):
        rf_sd.sort()

        rf = Reader_output("../outputs/" + rf_sd[ii]) 
        stdev = Reader_output("../outputs/" + rf_sd[len(rf_sd)//2 + ii]) 

        rf_dist = rf_sd[ii].split("_")
        stdev_dist = rf_sd[len(rf_sd)//2 + ii].split("_")

        if rf_dist[-2] == "OffDiag":
            title= f"Return Probability (right) and Standard Deviation (left) with a perturbation following a\n {rf_dist[1]}, W={rf[0][0]} and off the diagonal of the hamitonian."

        else:
            title= f"Return Probability (right) and Standard Deviation (left) with \na perturbation following a {rf_dist[1]} and W={rf[0][0]}."

        plot_2graphs(rf, stdev, [[0,1,2]], [[0,1,2]], 
                    [rf_sd[ii], rf_sd[len(rf_sd)//2 + ii]], axis=[["time", "return probability"], ["time", "standard deviation"]])

if grdst:    
    print("Ground State Amplitude in each lattice for differente perturbation values:")
    for ii in tqdm(grdst):
        grdst_dt = Reader_output("../outputs/"+ii)
        plot_graph_grdst(grdst_dt, [[0,1,2]], ii)
    
print("IPR value for differnt values of W in Aubry-Andre model:")
IPR_dt = Reader_output("../outputs/Wvs.IPR_QuasiPeriodicPert_withAlpha.txt")
print(IPR_dt[1])
lst = IPR_sorter(IPR_dt)
plot_graphs_IPR(lst, [[0,1,2,3]], "Wvs.IPR_QuasiPeriodicPert_withAlpha.txt")

if grdst_PBC:
    print("Ground State Amplitude in each lattice for differente perturbation values in PBC's:\n")
    for ii in tqdm(grdst_PBC):
        grdst_PBC_dt = Reader_output("../outputs/"+ii)
        print(ii)
        if len(grdst_PBC_dt)==3:
            plot_graph_grdst(grdst_PBC_dt, [[0]], style="-")
        else:
            plot_graph_grdst(grdst_PBC_dt, [[0, 1, 2]], style="-")

if S:    
    print("Entro:")
    for ii in tqdm(range(len(S))):
        S.sort()
        S_dt = Reader_output("../outputs/" + S[ii])
        S_dist = S[ii].split("_")
        plot_graph(S_dt, [[0,1,2]], S[ii], axis=["time", "entanglement entropy"])

if imbalance:
    for ii in tqdm(range(len(imbalance))):
        imbalance.sort()
        imb = Reader_output("../outputs/" + imbalance[ii])
        imb_dist = imbalance[ii].split("_")
        plot_graph(imb, [[0,1,2]], imbalance[ii],axis=["time", "imbalance"])


if IPR:
    print("IPR value evolution in time:")
    for ii in tqdm(range(len(IPR))):
        IPR.sort()
        IPR_dt = Reader_output("../outputs/" + IPR[ii])
        IPR_dist = IPR[ii].split("_")
        plot_graph(IPR_dt, [[0,1,2]], IPR[ii], axis=["time", "IPR"])


"""
for ii in tqdm(ESpec):
    ESpec_dt = Reader_output("../outputs/"+ii)
    plot_graph_ESpec(ESpec_dt, ii, style="o")


if grdst_PBC:
    print("Ground State Amplitude in each lattice for differente perturbation values in PBC's:\n")
    for ii in tqdm(grdst_PBC):
        grdst_PBC_dt = Reader_output("../outputs/"+ii)
        print(ii)
        if len(grdst_PBC_dt)==3:
            plot_graph_grdst(grdst_PBC_dt, [[0]], style="-")
        else:
            plot_graph_grdst(grdst_PBC_dt, [[0, 1, 2]], style="-")

if grdst_PBC:    
    print("Ground State Amplitude in each lattice for differente perturbation values in PBC's:\n")
    for ii in tqdm(grdst_PBC):
        grdst_PBC_dt = Reader_output("../outputs/"+ii)
        plot_graph_grdst(grdst_PBC_dt, [[0,1,2]], style="-")

if grdst:    
    print("Ground State Amplitude in each lattice for differente perturbation values:\n")
    for ii in tqdm(grdst):
        grdst_dt = Reader_output("../outputs/"+ii)
        plot_graph_grdst(grdst_dt, [[0,1,2]], ii)






if perialp:
    print("Return Probability and Standard Deviation for a periodic perturbation with a random phase deviation:\n")
    for ii in tqdm(range(len(perialp)//2)):
        perialp.sort()

        rf = Reader_output("../outputs/" + perialp[ii]) 
        stdev = Reader_output("../outputs/" + perialp[len(perialp)//2 + ii]) 

        rf_dist = perialp[ii].split("_")
        stdev_dist = perialp[len(perialp)//2 + ii].split("_")
        print(rf_dist, stdev_dist)

        if rf_dist[-2] == "OffDiag":
            title= f"Return Probability (right) and Standard Deviation (left) with a periodic perturbation\n, W={rf[0][0]}, a random phase and off the diagonal of the hamitonian."

        else:
            title= f"Return Probability (right) and Standard Deviation (left) with a periodic perturbation,\n W={rf[0][0]} and a random phase."

        plot_2graphs(rf, stdev, [[0,1,2]], [[0,1,2]], 
                    [perialp[ii], perialp[len(perialp)//2 + ii]], 
                    title=title, style="x")

if peri:
    print("Return Probability and Standard Deviation for a periodic perturbation:\n")
    for ii in tqdm(range(len(peri)//2)):
        peri.sort()

        rf = Reader_output("../outputs/" + peri[ii]) 
        stdev = Reader_output("../outputs/" + peri[len(peri)//2 + ii]) 

        rf_dist = peri[ii].split("_")
        stdev_dist = peri[len(peri)//2 + ii].split("_")
        print(rf_dist, stdev_dist)

        if rf_dist[-2] == "OffDiag":
            title= f"Return Probability (right) and Standard Deviation (left) with a periodic perturbation\n, W={rf[0][0]} and off the diagonal of the hamitonian."

        else:
            title= f"Return Probability (right) and Standard Deviation (left) with a periodic perturbation\n, and W={rf[0][0]}."

        plot_2graphs(rf, stdev, [[0,1,2]], [[0,1,2]], 
                    [peri[ii], peri[len(peri)//2 + ii]], 
                    title=title, style="x")


"""