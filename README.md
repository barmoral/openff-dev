Calculation and Analysis of Osmotic Values
================

Calculations can be performed with two different types of restraints: [Harmonic Potential (HP) and Flat-Bottom Potentials (FBPs)](osmotic_theory_summary.md).

### Simulation Details
* Python codes were written for simulation boxes with dimension of 4.8nm x 4.8nm x 14.4 nm.
* The number of ions in each simulation changes depending on the desired concentration in the box.
* The number of water molecules stays fixed at 11100.
* The sample simulations performed in this repository were equilibrated for 2 ns and their production run lasted 20 ns.

* Input coordinate files are generated using packmol. The packmol inputs can be found on the directory called ['structures'](structures/).

* The force constant (k, with units of kJ/mol/nm^2) used for HP and FBP restraints is different. While the FBP method is not very sensitive to the force constant value (in here I use the same value as the method developers, [Luo & Roux](https://pubs.acs.org/doi/10.1021/jz900079w)), HP is highly sensitive, and this value has to be estimated beforehand using the python notebook called ['force_constant_hps.ipynb'](https://github.com/barmoral/openff-dev/blob/main/structures/force_constant_hps.ipynb).

* I included example script files to use in a cluster called [run_HP.sh](https://github.com/barmoral/openff-dev/blob/main/run_HP.sh) and [run_FBP.sh](https://github.com/barmoral/openff-dev/blob/main/run_FBP.sh), which can also give an idea of how to run the calculations from any CLI.


### Steps to run Harmonic Potential calculations and analysis:
1. Define maximum concentration desired and corresponding number of ions needed in the given volume. (I used 4 molal as maximum concentration, which requires 217 NaCl molecules (217 of each ion)). 
2. Estimate the force constant using notebook [force_constant_hps.ipynb](https://github.com/barmoral/openff-dev/blob/main/force_constant_hps.ipynb). (0.68095403 kJ/mol/nm^2 was the result I obtained)
3. Run in cluster using script file [run_HP.sh](https://github.com/barmoral/openff-dev/blob/main/run_HP.sh). (I indicated the variables that need to be modified. Make sure the variable 'RESTRAINT' is set to HP)
4. The script file will run the calculation code, [osmotic_sim_dispatch.py](https://github.com/barmoral/openff-dev/blob/main/osmotic_sim_dispatch.py). If modifications to simulation set-up want to be done, they have to be done directly in the code.
5. The calculation code will generate three directories: 1) equilibration run using NVT ensemble input and output files for the simulation of the specified concentration and the corresponding replicate number (e.g., equil_sim_NVT_01m_r0) 2) equilibration run using anisotropic NPT ensemble input and output files for the simulation of the specified concentration and the corresponding replicate number (e.g., equil_sim_NPT_01m_r0) 3) production run using anisotropic NPT ensemble input and output files for the simulation of the specified concentration and the corresponding replicate number (e.g., prod_sim_01m_r0)
6. When all replicates for a concentration have finished running, the calculation code will generate a directory called 'result_files' for that specified concentration. In there it will copy the resulting topology file (.pdb) from the production run, and will convert its resulting trajectory .dcd file into a .xtc file and copy it there as well.  
7. Run the analysis notebook called [HP_analysis.ipynb](https://github.com/barmoral/openff-dev/blob/main/HP_final/HP_analysis.ipynb). Make sure to specify simulation set-up details and check that the directory names are appropiately formated at the beginning of the notebook (e.g., HP_MgCl_TIP3P/result_files_1m). (Note: In this code, a bootstrapping calculation will be a bottleneck, taking ~11 minutes to run if you use 500 bootstraps). 
8. Results and plots will automatically be saved to the 'result_files' directories holding the input files.


### Steps to run Flat-Bottom Potentials calculations and analysis:
1. Define which concentrations are wanted and their corresponding number of ions for the given volume. (I used 1 molal, 2 molal, and 3 molal, which respectively required 65, 128, and 188 NaCl molecules)
2. Define desired force constant. (I use 4184 kJ/mol/nm^2, which is the one [Luo & Roux](https://pubs.acs.org/doi/10.1021/jz900079w) recommended in their 2010 paper I was attempting to replicate)
3. Run in cluster using script file [run_FBP.sh](https://github.com/barmoral/openff-dev/blob/main/run_FBP.sh). (I indicated the variables that need to be modified. Make sure the variable 'RESTRAINT' is set to FBP)
4. The script file will run the calculation code, [osmotic_sim_dispatch.py](https://github.com/barmoral/openff-dev/blob/main/osmotic_sim_dispatch.py). If modifications to simulation set-up want to be done, they have to be done directly in the code.
5. The calculation code will generate three directories: 1) equilibration run using NVT ensemble input and output files for the simulation of the specified concentration and the corresponding replicate number (e.g., equil_sim_NVT_01m_r0) 2) equilibration run using anisotropic NPT ensemble input and output files for the simulation of the specified concentration and the corresponding replicate number (e.g., equil_sim_NPT_01m_r0) 3) production run using anisotropic NPT ensemble input and output files for the simulation of the specified concentration and the corresponding replicate number (e.g., prod_sim_01m_r0)
6. When all replicates for a concentration have finished running, the calculation code will generate a directory called 'result_files' for that specified concentration. In there it will copy the resulting topology file (.pdb) from the production run, and will convert its resulting trajectory .dcd file into a .xtc file and copy it there as well.  
7. Run the analysis notebook called [FBP_analysis.ipynb](https://github.com/barmoral/openff-dev/blob/main/FBP_final/FBP_analysis.ipynb). Make sure to specify simulation set-up details and check that the directory names are appropiately formated at the beginning of the notebook (e.g., FBP_NaCl_TIP3P/result_files_1m). 
8. Results and plots will automatically be saved to the 'result_files' directories holding the input files.


### Code to compare FBPs and HP results with experiments:
Run the notebook [results_FBP_vs_HP.ipynb](https://github.com/barmoral/openff-dev/blob/main/results_FBP_vs_HP.ipynb) which merges the results for both calculation types and plots them along with experiments to be able to compare both methods. A directory for the plots will be generated and plots will automatically save there. 