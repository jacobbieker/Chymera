1. hsct-poly:
  - Open hscf.in and change the parameters to what you want
  - compile hscf.poly with gfortran, unless there is already an a.out file in the directory
  - run ```./a.out < hscf.in``` to run the program
  - it should output files, what needs to be copied into the run directory is the fort.93 file, rename it to fort.2

2. Chymera point mass model with wiggle
  - Go to the ```bin``` folder, and open ```Makefile.DEFS``` 
  - Change the exe option to the run directory you need
  - Make sure Wiggle = 1 and Wiggle_restart = 0
  - run ```make clean``` then ```make```
  - change directory to the run directory
  - make sure the ```fort.5``` file, the ITYPE = 7
  - ```cd``` up to the directory with qsub, and run the job
  - After the initial setup, go back into ```Makefile.DEFS``` and change Wiggle_Restart = 1
  - Recompile with ```make clean``` then ```make```
  - Go into the run directory
  - ```cp starry_restart.TIMESTEP starry_restart.dat```
  - ```cp mfrac.dat mfrac.TIMESTEP```
  - ```cp saved.TIMESTEP fort.7```
  - In fort.5 change ITYPE = 1
  - run job
