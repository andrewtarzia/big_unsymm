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

from env_set import meta_path, calc_path


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
        smiles='NCCN',
        functional_groups=(stk.PrimaryAminoFactory(
            bonders=(1, ),
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

    logging.info(f'building metal')
    sqpl_unopt = stk.ConstructedMolecule(
        topology_graph=stk.metal_complex.CisProtectedSquarePlanar(
            metals=pd,
            ligands=bidentate,
        ),
    )
    sqpl_unopt.write(meta_unopt)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
