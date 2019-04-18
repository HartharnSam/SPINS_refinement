# SPINS grid refinement tools
This repository contains two main tools to make working with [SPINS](https://git.uwaterloo.ca/SPINS/SPINS_main) easier and more efficient: driver scripts to automate your workflow when sweeping parameter space, and a script to interpolate existing SPINS output to a higher resolution grid and restart the simulation. These tools are specifically written to use on the Compute Canada [Graham](https://docs.computecanada.ca/wiki/Graham) system. The SPINS grid refinement script should be quite portable between systems, but the driver scripts are especially tailored to Graham's scheduling software.

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
Grid refinement is done through the script `spins_refinement.m` in the `matlab` folder. This script has a lot in it, but it is easy to use. The typical steps you want to use to use the refinement script are:

 1. Adjust the parameters at the top of the script that specify the new grid size, which output to refine, and the interpolation method. The parameters are:

```matlab
 nx = 4;              % nx times of original resolution in x
 nz = 2;              % nz times of original resolution in z
 ii = 100;             % number of output to interpolate
 method = 'nearest';       % Matlab built-in interpolation method
                           % 'nearest', 'linear', 'spline' or 'cubic'
 ```
 The parameters `nx = 4` and `nz = 2` mean that the number of grid points in the horizontal and vertical directions will be increased by a factor of 4 and 2 respectively. The line ii = 100 determined what output would the user would like to refine and begin the new simulation from. This output will be the new initial condition for the higher resolution version of your case. Be aware that the output you choose must be sufficiently early that any transient motion due to the interpolation can be removed by viscosity or filtering before the motion of interest is analyzed. Lastly, you may also edit the variable ‘method’, which changes the interpolation method used. Currently, the only methods available are nearest, linear, spline, or cubic. In the future, a separate script for spectral interpolation may be added.
 2. Run the edited `spins_refinement.m` in the working directory. This script writes the interpolated fields to file in a new directory called `high_res`. It also writes a `spins.conf` file with the new grid specifications.
 3. Copy the submit script and `wavereader.x` into the `high_res` directory. You probably want to increase the number of processors and/or increase the allocated time in the submission script. For example, if you have `nx = 4` and `nz = 2`, you might increase number of processors 4 times, and double the allocated time. Now you can submit your high resolution simulation!

 If you want to change any parameters in `spins.conf` other than the resolution, you need to add an extra line to `spins_refinement.m`. For example, say you are interested in also changing the viscosity for the higher resolution case. Then to change the viscosity to 1e-7, you would add
 ```matlab
 params.visco = 1.e-7
 ```
 
 before the new `spins.conf` is written. A logical place could be after the new `Nx` and `Nz` are written:
 ```matlab

 % Update grids in spins.conf file
 params.Nx = Nx_new;
 params.Nz = Nz_new;
 params.visco = 1.e-7; % Add this
 ```

Note: `spins_refinment.m` requires the [SPINSmatlab](https://git.uwaterloo.ca/ddeepwel/SPINSmatlab) repository to run. Clone the SPINSmatlab repository and follow the instructions over there to add to your matlab path.