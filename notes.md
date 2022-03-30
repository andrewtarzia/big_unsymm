working directory: /data/atarzia/projects/big_unsymm
env_set.py sets environment, incl. directories for output

build_ligand.py
    build ligand precusor for cage construction -- includes optimisation
    generates ligands lowest energy conformer

build_metal_precusor.py
    will build the cis-protected Pd precursor for subsystems -- includes optimisation

build_cage.py
    will build 4 cage symmetries -- includes optimisation

build_subsystems.py
    will build triangles, squares and 6-mem rings from metal precursor and ligand -- includes optimisation

utilities.py
    defines utilities for construction
    same as unsymm basically.

optimisation.py
    defines optimisation sequence.
    same as unsymm.
    actually, more similar to cubism, with longer MD.

topologies.py
    contains the new topology graph for the ring and updated M12L24

analyse_cages.py
    for all desired systems:
        calculates total energies
        calculates metal-atom order parameter
        calculates ligand strain energy

analyse_subsystems.py
    for all desired systems:
        calculates total energies
        calculates metal-atom order parameter
        calculates ligand strain energy


plotting data...
