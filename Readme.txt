In order to run this code you only have to open the terminal inside this directory, 
have the eigen library installed in your computer and run the make file in the following 
way:

$ make TARGET=x PROGRAM=x

being x the name of the desired system to simulate. 

If the makefile is giving error because of the Eigen library, change the makefile line
-I/usr/include/eigen3 to -I/path/to/eigen3. You also need to have:
-/inc:
    functions.h
    iofile.h
-/src:
    functions.cpp
    iofile.cpp
    plotter.py

Here follows a list of all the programs possible to run:

x:
-> single_e-: Simulates a unidimensional lattice with an eletrons in the middle of it 
and a determined perturbation amplitude. The sizes of the lattice(L), the values of perturbation
amplitude(W), the nature of the hamiltonian (if the random field is in the diagonal or off 
diagonal) is asked as input at the begin of the program. You can insert more than one value of
L and W. This program calculates the evolution in time of the return probability of this
lattice point and the respective standard deviation.

-> aubry_andre_model: Simulates a unidimensional lattice with a periodic potential with an 
eletrons in the middle of it and for a determined perturbation amplitude. 
The sizes of the lattice(L), the values of perturbation amplitude(W), the nature of the 
hamiltonian (if the potential field is in the diagonal or off diagonal), if you want periodic 
boundary conditions is asked as input at the begin of the program is asked as input at the 
begin of the program. You can insert more than one value of L and W. This program calculates
the evolution in time of the return probability of this lattice point and the respective 
standard deviation.

-> aubry_andre_model_GS: Same system as the above program but it calculates the ground state 
in each lattice for different values of W and L.

-> Hofstadter_Butterfly: Same system as the one in the above program but it calculates the
energy spectrum of the system with different values of a  incommensurate parameter between 
0 and 1. This for a certain value of L and different values o W.

-> IPRvsW: Same system as the one in the above program but it calculates the IPR value 
denpendent of W for different L valeus specified by the user.

-> chain_e-: Simulates a unidimensional lattice with a random potential field in the 
hamiltonian diagonal with eletrons in the even indices of the lattice for a determined 
perturbation amplitude. This program calculates the temporal evolution of the systems 
imbalance and its entanglement entropy.

These programs output a .txt file in which there is all the needed information to plot These
quantities using the main.py program. To create the plots just move to the main directory 
and run the main.py program. The plot figures are saved in the plots directory.
