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

from env_set import meta_path, calc_path, gulp_path


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

    bidentate = stk.BuildingBlock(
        smiles='NCCCN',
        functional_groups=(stk.PrimaryAminoFactory(
            deleters=(),
        ), ),
    )
    pd = stk.BuildingBlock(
        smiles='[Pd+2]',
        functional_groups=(
            stk.SingleAtom(stk.Pd(0, charge=2))
            for i in range(4)
        ),
        position_matrix=[[0, 0, 0]],
    )

    meta_unopt = os.path.join(_wd, 'meta_unopt.mol')
    meta_opt = os.path.join(_wd, 'meta_opt.mol')

    logging.info(f'building metal')
    sqpl_unopt = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
            metals=pd,
            ligands=bidentate,
            optimizer=stk.MCHammer(
                num_steps=500,
                target_bond_length=1.8,
            ),
        ),
    )
    sqpl_unopt.write(meta_unopt)

    if not os.path.exists(meta_opt):
        CG = False
        logging.info(f'UFF4MOF optimisation 1 of meta CG: {CG}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=gulp_path(),
            maxcyc=1000,
            metal_FF={46: 'Pd4+2'},
            metal_ligand_bond_order='',
            output_dir=os.path.join(_cd, f'meta_gulp'),
            conjugate_gradient=CG,
        )
        gulp_opt.assign_FF(sqpl_unopt)
        sqpl_gulp = gulp_opt.optimize(mol=sqpl_unopt)
        sqpl_gulp.write(meta_opt)
    else:
        sqpl_gulp = sqpl_unopt.with_structure_from_file(meta_opt)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
