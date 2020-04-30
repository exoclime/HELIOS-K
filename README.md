# HELIOS-K #
HELIOS-K is an opacity calculator, running on GPUs.
#### Authors: Simon Grimm, Kevin Heng ####

# Updates #
## version 1.67 ##
The `Molecule` and `useHITEMP` arguments in the `param.dat` file are not valid anymore. Species must now be set by `Species Name`. The database is written from the `< species >.param` file.

## version 1.65 ##
All species must have now a `< species >.param` file, which contains all necessary information about the line list.

