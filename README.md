# RPCOVA (A relief from unmeasured confounding, no instrument or negative control)
Causal inference from nonrandomized studies using Resistant Population Calibration Of VAriance


The template code to implement the method is in rpcova.R. Please run the code line by line. This initial commit does not provide a function.

The constrained optimazation is done using the Gurobi solver. See direction for installation in https://www.gurobi.com/

Results in Sections 5.2 and 5.3 require longer computation times. The code provided here is set up to run in a cluster as a slurm job. Please run the .sbatch files within each folder to generate the results files. These results are then compiled to create the figures and tables.

Please contact us if you have difficulty using the code. Email: Bikram Karmakar (bkarmakar@ufl.edu) Zikun Qin (qinzk96@hotmail.com and qinzikun@ufl.edu)
