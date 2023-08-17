# NUCLEOSOME_SEQUENCE_DEPENDENCE
Fortran 90 programs for analyzing various stuff of nucleosome. Most of the programs will be used to directly analyse trajectories generated from MD engine. However, these programs need dcd formatted trajectories along with a reference pdb file that will be basically needed for selection. Rest of the programs will be used to post-process Curves+ (an open source program for analysizing nucleosome deformation variables such as Rise, Shift, Slide, Twist, Roll, Tilt) generated data for calculating deformation energy cost and deformation force constant. The sample input files for running each of the programs are in INPUT folder with the same name as that of the programs with an extension .inp. 

Compiler requirement: gfortran/f95
Usage: 
compilation:   f95/gfortran program_name.F90
Execution:     ./a.out<program_name.inp>program_name.log &

Specific details of all the programs and the requirements are mentioned in the preambles of each of the programs


Note: The inputs provided in the INPUT folder are samples only. Users need to modify the input file as per their need. If there is any issue to run the programs, please contact me at prabir07chem@gmail.com
