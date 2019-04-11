# SPINS grid refinement tools

## Refinement script

## Tutorial

## Workflow


Thoughts about writing tutorial:

  * Add grid refinement tutorial to the end of the existing SPINS tutorial.
  * Use Will's code to sweep parameter space
  * Link to grid refinement code on git repo
  * Workflow diagram + some instuctions on how to use the grid refinement code here.

`spins_refinement.m` workflow:


Modifying parameters:
 * Modify tilttank.txt, driver_tilttank.txt to control Lx, Lz, Nx, Nz, density difference and "amplitude": how much fluid is behind the gate.
 * Go into `matlab.sh`, change casenames string to match those in text file
 * Run matlab.sh