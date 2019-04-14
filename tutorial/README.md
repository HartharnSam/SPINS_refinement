# SPINS tutorial extension

This tutorial is an extension of the [SPINS tutorial](https://wiki.math.uwaterloo.ca/fluidswiki/index.php?title=SPINS_Tutorial) on the fluids wiki. In that tutorial, you manually changed parameters in the spins.conf file to run different cases. In this tutorial you will:

1. use a text file specifying parameters for a selection of cases to automatically set up the cases and submit them all to the compute canada scheduler.
2. use the matlab script `spins_refinement.m` to refine the grid for one of the above cases and restart the simulation at higher resolution from a specified timestep.

Note that this extension uses wavereader instead of modifying a `.cpp` case file.

These tools are intended to help us in carrying out well designed suites of simulations. The idea is that once you have an idea for an experiment in mind, the general workflow will be to:

* Design a base case to see if what you want to do is possible. Run this at low resolution as a proof of concept.
* Once you have the base case running, design a set of cases to explore the effect of certain parameters, or find the ideal parameter value for your experiment. These tools automate the process of creating the case files, directories, and submitting the jobs.
* Restart one or more of the above cases from a specified time in the simulation at higher resolution. For example, let the wave propagate and form at low resolution, and increase the resolution before it begins shoaling.

## Setup instructions
Clone the [SPINS refinement git repository](https://git.uwaterloo.ca/tghill/SPINS-refinement/) and copy the `tutorial` directory to the directory you have been running SPINS from.

## Designing cases to explore parameter space

### Using `mode1_mode2.txt`
Some of the parameters that might be interesting to change are `H1`, `H2`, `L1`, `L2`, and `pyc_loc`.

The text file `mode1_mode2.txt` controls the cases we will run. Open this file and look at its contents. There are 5 cases, each with a different set of parameters:

1. base. This is the base case, with the default parameters as in the original spins.conf file for this tutorial.
2. mode1. This case changes the mode 2 amplitude and length to 0, meaning we will only have a mode 1 wave propagating to the right.
3. mode2. This case changes the mode 1 amplitude and length to 0, meaning we will only have a mode 2 wave propagating to the right.
4. gravity. This cases uses the parameters as described above to create a gravity current.
5. lamp. This cases increases the height and width of the mode 1 and mode 2 regions to create larger amplitude waves. The pycnocline is also set to be in the middle of the tank.

The matlab script `driver_mode1_mode2.m` sets up the cases. It creates directories for each case (named according to the casename), writes the parameters to a new `spins.conf` file, and creates the initial u, v, w, and rho files. This script has its own documentation on the git repository, and the code is well commented.

Now we actually run these cases.

First, open matlab and run `driver_mode1_mode2.m` to set up the cases. For now, look in each directory it created so you can see what it put in each directory.

Now, submit the `matlab.sh` script to the scheduler (`sbatch matlab.sh`). This job runs `wavereader.x` for each case, and submits a job for each case. Specifically, it submits the `submit.sh` script from each case directory.


### What to do to make your own example
Or should this be in the package documentation on git? Probably there.


## Increasing resolution


## To do still:

 * make directory `examples` with a selection of tables and corresponding driver files
 * Test the instructions here and make sure everything runs.
 * Let Andrew read it over (Monday?) then put it on the wiki page instead of the git page.
