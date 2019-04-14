# SPINS tutorial extension

This tutorial is an extension of the [SPINS tutorial](https://wiki.math.uwaterloo.ca/fluidswiki/index.php?title=SPINS_Tutorial) on the fluids wiki. In that tutorial, you manually changed parameters in the spins.conf file to run different cases. In this tutorial you will:

1. use a text file specifying parameters for a selection of cases to automatically set up the cases and submit them all to the compute canada scheduler.
2. use the matlab script `spins_refinement.m` to refine the grid for one of the above cases and restart the simulation at higher resolution from a specified timestep.

These tools are intended to help us in carrying out well designed suites of simulations. The idea is that once you have an idea for an experiment in mind, the general workflow will be to:

* Design a base case to see if what you want to do is possible. Run this at low resolution as a proof of concept.
* Once you have the base case running, design a set of cases to explore the effect of certain parameters, or find the ideal parameter value for your experiment. These tools automate the process of creating the case files, directories, and submitting the jobs.
* Restart one or more of the above cases from a specified time in the simulation at higher resolution. For example, let the wave propagate and form at low resolution, and increase the resolution before it begins shoaling.

## Designing cases to explore parameter space
