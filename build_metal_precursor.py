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

from env_set import (
    meta_path,
    calc_path,
    xtb_path,
    gulp_path,
)
from utilities import AromaticCNCFactory
from spinner import rotate_fgs


def get_chosen_conformer(mol):
    """
    Select and optimize a conformer with desired directionality.

    Currently:
        Best directionality will be defined by the smallest
        N-ligand centroid-N angle.

    """

    lowe_mol = stk.BuildingBlock.init_from_molecule(
        molecule=mol,
        functional_groups=(stk.BromoFactory(), ),
    )

    out_molecule = rotate_fgs(lowe_mol)
    return out_molecule


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

    if not os.path.exists(_cd):
        os.mkdir(_cd)

    pyridine = stk.BuildingBlock(
        smiles='C1=CC=NC=C1',
        functional_groups=(AromaticCNCFactory(), ),
    )
    meta = stk.BuildingBlock(
        smiles='C1=CC(=CN=C1)Br',
        functional_groups=(AromaticCNCFactory(), ),
    )
    para = stk.BuildingBlock(
        smiles='C1=CN=CC=C1Br',
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

    _topos = {
        'cis1m1p': {
            pyridine: (0, 1),
            para: (2, ),
            meta: (3, ),
        },
        'cis2m0p': {
            pyridine: (0, 1),
            meta: (2, 3),
        },
        'cis0m2p': {
            pyridine: (0, 1),
            para: (2, 3),
        },
        'trans1m1p': {
            pyridine: (0, 2),
            para: (1, ),
            meta: (3, ),
        },
        'trans2m0p': {
            pyridine: (0, 2),
            meta: (1, 3),
        },
        'trans0m2p': {
            pyridine: (0, 2),
            para: (1, 3),
        },
    }

    for topo in _topos:
        meta_unopt = os.path.join(_wd, f'{topo}_unopt.mol')
        meta_rot = os.path.join(_wd, f'{topo}_rot.mol')
        meta_gulp = os.path.join(_wd, f'{topo}_gulp.mol')
        meta_opt = os.path.join(_wd, f'{topo}_opt.mol')

        sqpl_unopt = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.SquarePlanar(
                metals=pd,
                ligands=_topos[topo],
                optimizer=stk.MCHammer(
                    num_steps=500,
                    target_bond_length=1.8,
                ),
            ),
        )
        sqpl_unopt.write(meta_unopt)

        if '2m' in topo:
            if os.path.exists(meta_rot):
                sqpl_rot = sqpl_unopt.with_structure_from_file(
                    path=meta_rot,
                )
            else:
                logging.info(f'getting optimal conformer of {topo}')
                sqpl_rot = get_chosen_conformer(mol=sqpl_unopt)
                sqpl_rot.write(meta_rot)
        else:
            sqpl_rot = stk.BuildingBlock.init_from_molecule(sqpl_unopt)

        # Do Gulp opt first.
        if not os.path.exists(meta_gulp):
            CG = False
            logging.info(f'UFF4MOF optimisation 1 of {topo} CG: {CG}')
            gulp_opt = stko.GulpUFFOptimizer(
                gulp_path=gulp_path(),
                maxcyc=1000,
                metal_FF={46: 'Pd4+2'},
                metal_ligand_bond_order='',
                output_dir=os.path.join(_cd, f'{topo}_gulp'),
                conjugate_gradient=CG,
            )
            gulp_opt.assign_FF(sqpl_rot)
            sqpl_gulp = gulp_opt.optimize(mol=sqpl_rot)
            sqpl_gulp.write(meta_gulp)
        else:
            sqpl_gulp = sqpl_rot.with_structure_from_file(meta_gulp)

        if not os.path.exists(meta_opt):
            logging.info(f'xtb optimisation of {topo}')
            xtb_opt = stko.XTB(
                xtb_path=xtb_path(),
                output_dir=os.path.join(_cd, f'{topo}_xtbopt'),
                gfn_version=2,
                num_cores=6,
                opt_level='normal',
                charge=2,
                num_unpaired_electrons=0,
                max_runs=1,
                calculate_hessian=False,
                unlimited_memory=True,
                solvent=None,
            )
            sqpl_opt = xtb_opt.optimize(mol=sqpl_gulp)
            sqpl_opt.write(meta_opt)
        else:
            sqpl_opt = sqpl_gulp.with_structure_from_file(meta_opt)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
