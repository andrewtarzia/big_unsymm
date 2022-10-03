# big_unsymm
Construct and study large (M12L24) cages with unsymmetrical ligands

DOI for publication: 10.1039/d2sc03856k

Installation:

Clone the repo and use enrivonment.yml and Anaconda or Miniconda to install the Python 3.9 environment including RDKIT and XTB from conda-forge.
CREST must be installed by the user (version 2.11.2 was used). In env_set.py, I set environment variables, including directories for output. E.g. the working directory for me is `/data/atarzia/projects/big_unsymm`.

Usage:

`build_ligand.py`: 
    build ligand precusor for cage construction -- includes optimisation. Generates ligands lowest energy conformer

`build_metal_precusor.py`:
    will build the cis-protected Pd precursor for subsystems -- includes optimisation. Not used in the manuscript.

`build_cage.py`:
    will build 4 cage symmetries -- includes optimisation

`build_subsystems.py`:
    will build triangles, squares and 6-mem rings from metal precursor and ligand -- includes optimisation. Not used in the manuscript.

`utilities.py`:
    defines utilities for construction

`optimisation.py`:
    defines optimisation sequences.

`topologies.py`:
    contains the new topology graph for the ring and updated M12L24
    The metal complex topology was not used in this work.

`analyse_cages.py`:
    for all desired systems:
    
        * calculates total energies
        * calculates metal-atom order parameter
        * calculates ligand strain energy

`analyse_subsystems.py` and `analyse_manual_subsytems.py`:
    for all desired systems: Not used in the manuscript.
    
        * calculates total energies
        * calculates metal-atom order parameter
        * calculates ligand strain energy
    


`plotting.py`:
    Utilities for plotting.

`spinner.py`:
    Utilities for using SpinDry to orient building blocks.
    Not used in the manuscript.

`m24_optimising.py`:
    Script for running the optimisation of Pd24L48 systems from XYZ file.
