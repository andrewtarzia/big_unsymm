#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build the metal precursor in this project.

Author: Andrew Tarzia

"""

import logging
import sys
import os
import stk
import stko

from env_set import calc_path, meta_path, gulp_path
from spinner import rotate_fgs
from utilities import AromaticCNCFactory
from topologies import OneTwoSquarePlanar


def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 0 arguments:'
        )
        sys.exit()
    else:
        pass

    _wd = meta_path()
    _cd = calc_path()

    if not os.path.exists(_wd):
        os.mkdir(_wd)

    bidentate = stk.BuildingBlock(
        smiles='C1=C(CCCC2=CC=CC=N2)N=CC=C1',
        functional_groups=(AromaticCNCFactory(), ),
    )
    monodentate1 = stk.BuildingBlock(
        smiles='C1=CN=CC=C1Br',
        functional_groups=(AromaticCNCFactory(), ),
    )
    monodentate2 = stk.BuildingBlock(
        smiles='C1=CC(=CN=C1)C2=CC=C(C=C2)Br',
        functional_groups=(AromaticCNCFactory(), ),
    )
    monodentate3 = stk.BuildingBlock(
        smiles='C1=CC=NC=C1',
        functional_groups=(AromaticCNCFactory(), ),
    )
    pd = stk.BuildingBlock(
        smiles='[Pd+2]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2))
            for i in range(4)
        ),
        position_matrix=[[0, 0, 0]],
    )

    # Define series of topologies to build.
    _topos = {
        'm1': {
            'tg': OneTwoSquarePlanar(
                metals=pd,
                ligands={
                    monodentate1: (0, ),
                    monodentate2: (1, ),
                    bidentate: (2, ),
                },
            ),
            'charge': 2*1,
        },
        'm2': {
            'tg': OneTwoSquarePlanar(
                metals=pd,
                ligands={
                    monodentate2: (0, 1),
                    bidentate: (2, ),
                },
            ),
            'charge': 2*1,
        },
        'm3': {
            'tg': OneTwoSquarePlanar(
                metals=pd,
                ligands={
                    monodentate1: (0, 1),
                    bidentate: (2, ),
                },
            ),
            'charge': 2*1,
        },
        'm4': {
            'tg': OneTwoSquarePlanar(
                metals=pd,
                ligands={
                    monodentate1: (1, ),
                    monodentate2: (0, ),
                    bidentate: (2, ),
                },
            ),
            'charge': 2*1,
        },
        't1': {
            'tg': stk.metal_complex.SquarePlanar(
                metals=pd,
                ligands={
                    monodentate1: (2, ),
                    monodentate2: (0, ),
                    monodentate3: (1, 3),
                },
            ),
            'charge': 2*1,
        },
        't2': {
            'tg': stk.metal_complex.SquarePlanar(
                metals=pd,
                ligands={
                    monodentate1: (0, 2),
                    monodentate3: (1, 3),
                },
            ),
            'charge': 2*1,
        },
        't3': {
            'tg': stk.metal_complex.SquarePlanar(
                metals=pd,
                ligands={
                    monodentate2: (0, 2),
                    monodentate3: (1, 3),
                },
            ),
            'charge': 2*1,
        },
        't4': {
            'tg': stk.metal_complex.SquarePlanar(
                metals=pd,
                ligands={
                    monodentate1: (0, ),
                    monodentate2: (2, ),
                    monodentate3: (1, 3),
                },
            ),
            'charge': 2*1,
        },
    }
    # Build them all.
    for topo in _topos:
        unopt_file = os.path.join(_wd, f'{topo}_unopt.mol')
        rot_file = os.path.join(_wd, f'{topo}_rot.mol')
        opt_file = os.path.join(_wd, f'{topo}_opt.mol')

        tg = _topos[topo]['tg']
        charge = _topos[topo]['charge']
        logging.info(f'building {topo}')
        unopt_mol = stk.ConstructedMolecule(tg)
        unopt_mol.write(unopt_file)

        # Do some forced ligand rotations.
        if 't' in topo:
            rot_mol = rotate_fgs(stk.BuildingBlock.init_from_molecule(
                molecule=unopt_mol,
                functional_groups=(stk.BromoFactory(), ),
            ))
            rot_mol.write(rot_file)
        else:
            rot_mol = stk.BuildingBlock.init_from_molecule(unopt_mol)

        if not os.path.exists(opt_file):
            logging.info(f'Gulp opt of {topo}')
            output_dir = os.path.join(_cd, f'{topo}_gulp')
            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=gulp_path(),
                maxcyc=300,
                metal_FF={46: 'Pd4+2'},
                metal_ligand_bond_order='',
                output_dir=output_dir,
                conjugate_gradient=True,
            )
            gulp_opt.assign_FF(rot_mol)
            gulp_mol = gulp_opt.optimize(mol=rot_mol)

            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=gulp_path(),
                maxcyc=300,
                metal_FF={46: 'Pd4+2'},
                metal_ligand_bond_order='',
                output_dir=output_dir,
                conjugate_gradient=False,
            )
            gulp_opt.assign_FF(gulp_mol)
            gulp_mol = gulp_opt.optimize(mol=gulp_mol)
            gulp_mol.write(opt_file)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
