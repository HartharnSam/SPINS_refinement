# SPINS grid refinement tools

## MATLAB driver scripts
The set of driver scripts (`driver_mode1_mode2.m`, `mode1_mode2.txt`, `matlab.sh`, and `submit.sh`) automates the workflow of submitting multiple jobs to explore parameter space. See the tutorial on the fluids wiki, using the files in the `tutorial/` folder.

The following section describes how to edit these scripts to work for your specific experiment.

### Writing the text file
First, decide what parameters will be case dependent, and which will be case independent. Write a text file (see format in `tutorial/mode1_mode2.txt`) specifying the parameters for each case you want to run.

### Writing the driver script
The driver script (see `tutorial/driver_mode1_mode2.m`) reads the parameters from the text file, writes a corresponding `spins.conf` file, and creates and writes the initial velocity and density fields. The case independent parameters should be specified at the top of the script as constants, and any case dependent parameters should be read in from the text table. You may have to change what is written to spins.conf

### Submitting the jobs to Graham
The `matlab.sh` submission script submits the jobs to the Graham scheduler. You must change the `casenames=...` line in this script to match the casenames specified in the text cases file. This script then submits a job from each case directory using `submit.sh`, which runs `wavereader.x`.

## SPINS grid refinement
Grid refinement is done through the script `spins_refinement.m` in the `matlab` folder.
