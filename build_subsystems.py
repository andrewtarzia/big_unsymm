#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build the subsystems in this project.

Author: Andrew Tarzia

"""

import logging
import sys
import os
import stk

from env_set import cage_path, calc_path, meta_path, liga_path
from utilities import AromaticCNCFactory
from optimisation import optimisation_sequence


def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 0 arguments:'
        )
        sys.exit()
    else:
        pass

    li_path = liga_path()
    lig_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(li_path, 'lig_opt.mol'),
        functional_groups=(AromaticCNCFactory(), ),
    )
    me_path = meta_path()
    corner_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 'meta_opt.mol'),
        functional_groups=(
            stk.SmartsFunctionalGroupFactory(
                smarts='[#46]~[#7]#[#6]~[#1]',
                bonders=(0, ),
                deleters=(1, 2, 3),
            ),
        ),
    )
    _wd = cage_path()
    _cd = calc_path()

    if not os.path.exists(_wd):
        os.mkdir(_wd)

    if not os.path.exists(_cd):
        os.mkdir(_cd)

    # Define series of topologies to build.
    _init_opts = stk.MCHammer(target_bond_length=1.5, num_steps=500)

    _topos = {
        # Triangles.
        'tri1': {},
        # Squares.
        'sqr1': {},
        # Hexagons.
        'hex1': {},
    }
    # Build them all.
    for topo in _topos:
        unopt_file = os.path.join(_wd, f'{topo}_unopt.mol')
        opt_file = os.path.join(_wd, f'{topo}_opt.mol')

        tg = _topos[topo]['tg']
        charge = _topos[topo]['charge']
        logging.info(f'building {topo}')
        unopt_mol = stk.ConstructedMolecule(tg)
        unopt_mol.write(unopt_file)
        continue
        if not os.path.exists(opt_file):
            logging.info(f'optimising {topo}')
            opt_mol = optimisation_sequence(unopt_mol, charge=charge)
            opt_mol.write(opt_file)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
