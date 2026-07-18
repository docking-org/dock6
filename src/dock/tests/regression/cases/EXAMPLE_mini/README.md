# EXAMPLE_mini — case template (NOT runnable as-is)

This is a layout template, not a working case. To make it real:

1. Put a seeded DOCK_GA input in `dock.in` (fixed `random_seed`, minimizer ON, few
   generations, small ensemble).
2. Put everything `dock.in` references into `inputs/` (receptor grid, fragment libraries,
   starting ligand mol2, parameter files).
3. List the scientific output basenames to compare in `manifest.txt`.
4. Capture goldens from the UNMODIFIED build:  ../../capture_golden.sh . /path/to/dock6
5. Rename this directory to something descriptive (e.g. `tournament_smallmol`).
