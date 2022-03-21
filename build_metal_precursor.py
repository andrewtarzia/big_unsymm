#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build the ligand in this project.

Author: Andrew Tarzia

"""

import logging
import sys
import os
import stk
import stko

from env_set import meta_path, calc_path, xtb_path
from utilities import AromaticCNCFactory


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

    meta_unopt = os.path.join(_wd, f'meta_unopt.mol')
    meta_opt = os.path.join(_wd, f'meta_opt.mol')

    pyridine = stk.BuildingBlock(
        smiles='C1=CC=NC=C1',
        functional_groups=(AromaticCNCFactory(), ),
    )
    amine = stk.BuildingBlock(
        smiles='C#N',
        functional_groups=(
            stk.SmartsFunctionalGroupFactory(
                smarts='[#6]~[#7]',
                bonders=(1, ),
                deleters=(),
            ),
        ),
    )
    pd = stk.BuildingBlock(
        smiles='[Pd+2]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2))
            for i in range(4)
        ),
        position_matrix=[[0, 0, 0]],
    )

    sqpl_unopt = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.SquarePlanar(
            metals=pd,
            ligands={amine: (0, 1),pyridine: (2, 3)},
            optimizer=stk.MCHammer(
                num_steps=300,
                target_bond_length=1.8,
            ),
        ),
    )
    sqpl_unopt.write(meta_unopt)

    if not os.path.exists(meta_opt):
        logging.info(
            f'xtb optimisation in {os.path.join(_cd, "meta_xtbopt")}'
        )
        xtb_opt = stko.XTB(
            xtb_path=xtb_path(),
            output_dir=os.path.join(_cd, 'meta_xtbopt'),
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
        sqpl_opt = xtb_opt.optimize(mol=sqpl_unopt)
        sqpl_opt.write(meta_opt)



if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
