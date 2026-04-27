### DOCK 6.13.1

You can access the tar/zip files of the source code in the [Releases](https://github.com/docking-org/dock6/releases) section at right -->

## Quick Overview
DOCK is molecular modeling program used to identify potential binding geometries and interactions of a molecule to a target. Specifically, docking is the identification of the low-energy binding modes of a small molecule, or ligand, within the active site of a macromolecule, or receptor, whose structure is known. A compound that interacts strongly with, or binds, a receptor associated with a disease may inhibit its function and thus act as a drug. Solving the docking problem computationally requires an accurate representation of the molecular energetics as well as an efficient algorithm to search the potential binding modes.

Historically, the DOCK algorithm addressed rigid body docking using a geometric matching algorithm to superimpose the ligand onto a negative image of the binding pocket. Important features that improved the algorithm's ability to find the lowest-energy binding mode, including force-field based scoring, on-the-fly optimization, an improved matching algorithm for rigid body docking and an algorithm for flexible ligand docking, have been added over the years. For more information on past versions of DOCK, click here.

With the release of DOCK 6, we continue to improve the algorithm's ability to predict ligand binding poses by adding new features like force-field scoring, enhanced solvation models, reference-based scoring options, and de novo design. For more information about the current release of DOCK, click here.


## What can DOCK6 do for you?
We and others have used DOCK for the following applications:

- predict binding modes of small molecule-protein complexes
- search databases of ligands for compounds that mimic the inhibitory binding interactions of an experimentally validated inhibitor
- search databases of ligands for compounds that bind a particular site of a specific protein
- search databases of ligands for compounds that bind nucleic acid targets
- examine possible binding orientations of protein-protein and protein-DNA complexes
- help guide synthetic efforts by examining small molecules that are computationally derived
- many more...

## New to DOCK6.13:

New methodologies have been added to DOCK_DN that allow for users to bias the selection of both fragments and torsions toward those of higher frequency in the provided set. This provides users with finer control over the fragment and torsion compositions of their final ensembles. This method can be enabled in a standard DOCK_DN run with no additional processing or input files, so long as the libraries provided have associated frequencies that would be output with standard fragment library generation in DOCK6. [Bickel et al., J.Comput. Chem. 2024](https://onlinelibrary.wiley.com/doi/10.1002/jcc.27508)

Filtering molecules in DOCK_GA by a soft molecular weight cutoff has been added. Now, users can allow a chance for molecules to be accepted beyond this cutoff, enabling some deviation around the cutoffs. Changes to mutation selection and how fragments are chosen for mutation have been modified such that DOCK_GA will no longer select fragments incompatible with the attempted mutation type.

Users can now toggle the use of Ligand Efficiency (introduced in DOCK6.12) when using Grid Score via Descriptor Score. This new feature can be used alongside Grid Score and will be calculated as: Ligand Efficiency = (Grid Score)/(# Active Heavy Atoms)

Users now have explicit control over the number of conformers DOCK stores during anchor-and-grow docking prior to clustering. The old parameter num_scored_conformers has been removed and replaced with two new parameters num_final_scored_poses and num_preclustered_conformers. Previously, this control was grouped under a single parameter that also controlled the final number of scored poses written out.

## Documentation and Installation
The documentation is online at `http://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm` and located in this repository in `~/docs/dock6_manual.html`

For complete installation instructions, including the use of other compilers, RDKit, and Docker, see the installation section of the User's Manual.

To start, you can obtain this code by cloning the repository (or getting a previous release version from the [Releases page](https://github.com/docking-org/dock6/releases)).

`git clone https://github.com/docking-org/dock6.git`

Once unpacked, installation starts with `./configure -help` in the `~/install` directory.  That is, in a terminal window type this (and then press the "enter" key):

`cd install; ./configure -help`

For more information, including tutorials, bug fixes, etc., please consult
the UCSF DOCK Web page:

http://dock.compbio.ucsf.edu/


