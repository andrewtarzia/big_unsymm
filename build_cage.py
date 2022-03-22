#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all cages in this project.

Author: Andrew Tarzia

"""

import logging
import sys
import os
import stk

from utilities import AromaticCNCFactory
from env_set import cage_path, calc_path, liga_path
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
        path=os.path.join(li_path, 'lig_lowe.mol'),
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

    _wd = cage_path()
    _cd = calc_path()

    if not os.path.exists(_wd):
        os.mkdir(_wd)

    if not os.path.exists(_cd):
        os.mkdir(_cd)

    # Define series of topologies to build.
    # _init_opts = stk.MCHammer(target_bond_length=2.5, num_steps=1000)
    _init_opts = stk.NullOptimizer()
    _react_factory = stk.DativeReactionFactory(
        stk.GenericReactionFactory(
            bond_orders={
                frozenset({
                    stk.GenericFunctionalGroup,
                    stk.SingleAtom,
                }): 9,
            },
        ),
    )

    _topos = {
        'def': {
            'tg': stk.cage.M12L24(
                building_blocks={
                    pd: range(0, 12),
                    lig_bb: range(12, 36),
                },
                optimizer=_init_opts,
                reaction_factory=_react_factory,
            ),
            'charge': 2*12,
        },
        'oh': {
            'tg': stk.cage.M12L24(
                building_blocks={
                    pd: range(0, 12),
                    lig_bb: range(12, 36),
                },
                # Based on A-1.
                vertex_alignments={
                    # Right.
                    12: 0, 13: 0, 14: 1, 15: 1,
                    # Left.
                    16: 1, 17: 1, 18: 0, 19: 0,
                    # Top.
                    20: 0, 21: 0, 22: 1, 23: 1,
                    # Bottom.
                    24: 1, 25: 1, 26: 0, 27: 0,
                    # Front.
                    28: 0, 29: 1, 30: 0, 31: 1,
                    # Back.
                    32: 0, 33: 1, 34: 0, 35: 1,
                },
                optimizer=_init_opts,
                reaction_factory=_react_factory,
            ),
            'charge': 2*12,
        },
        's6': {
            'tg': stk.cage.M12L24(
                building_blocks={
                    pd: range(0, 12),
                    lig_bb: range(12, 36),
                },
                # Based on A-2.
                vertex_alignments={
                    # Right.
                    12: 0, 13: 1, 14: 0, 15: 1,
                    # Left.
                    16: 1, 17: 0, 18: 1, 19: 0,
                    # Top.
                    20: 0, 21: 0, 22: 1, 23: 1,
                    # Bottom.
                    24: 1, 25: 1, 26: 0, 27: 0,
                    # Front.
                    28: 0, 29: 1, 30: 0, 31: 0,
                    # Back.
                    32: 0, 33: 1, 34: 1, 35: 1,
                },
                optimizer=_init_opts,
                reaction_factory=_react_factory,
            ),
            'charge': 2*12,
        },
        'c4h': {
            'tg': stk.cage.M12L24(
                building_blocks={
                    pd: range(0, 12),
                    lig_bb: range(12, 36),
                },
                # Based on A-3.
                vertex_alignments={
                    # Right.
                    12: 1, 13: 0, 14: 1, 15: 0,
                    # Left.
                    16: 0, 17: 1, 18: 0, 19: 1,
                    # Top.
                    20: 0, 21: 1, 22: 0, 23: 1,
                    # Bottom.
                    24: 1, 25: 0, 26: 1, 27: 0,
                    # Front.
                    28: 0, 29: 0, 30: 1, 31: 1,
                    # Back.
                    32: 1, 33: 1, 34: 0, 35: 0,
                },
                optimizer=_init_opts,
                reaction_factory=_react_factory,
            ),
            'charge': 2*12,
        },
        'd2h': {
            'tg': stk.cage.M12L24(
                building_blocks={
                    pd: range(0, 12),
                    lig_bb: range(12, 36),
                },
                # Based on C-1.
                vertex_alignments={
                    # Right.
                    12: 0, 13: 0, 14: 1, 15: 1,
                    # Left.
                    16: 1, 17: 1, 18: 0, 19: 0,
                    # Top.
                    20: 0, 21: 0, 22: 1, 23: 1,
                    # Bottom.
                    24: 1, 25: 1, 26: 0, 27: 0,
                    # Front.
                    28: 1, 29: 0, 30: 1, 31: 0,
                    # Back.
                    32: 1, 33: 0, 34: 1, 35: 0,
                },
                optimizer=_init_opts,
                reaction_factory=_react_factory,
            ),
            'charge': 2*12,
        },
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

        if not os.path.exists(opt_file):
            logging.info(f'optimising {topo}')
            opt_mol = optimisation_sequence(
                mol=unopt_mol,
                name=topo,
                charge=charge,
                calc_dir=_cd,
            )
            opt_mol.write(opt_file)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
