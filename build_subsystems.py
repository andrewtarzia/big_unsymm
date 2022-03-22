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
from optimisation import optimisation_sequence
from m6l6 import M6L6
from utilities import AromaticCNC, AromaticCNCFactory


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
    me_path = meta_path()
    corner_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 'meta_unopt.mol'),
        functional_groups=(
            stk.SmartsFunctionalGroupFactory(
                smarts='[#46]~[#7]',
                bonders=(0, ),
                deleters=(),
                placers=(0, 1),
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
    _init_opts = stk.NullOptimizer()
    # _init_opts = stk.MCHammer(target_bond_length=2.5, num_steps=1000)
    _react_factory = stk.DativeReactionFactory(
        stk.GenericReactionFactory(
            bond_orders={
                frozenset({
                    stk.GenericFunctionalGroup,
                    AromaticCNC,
                }): 9,
            },
        ),
    )

    _topos = {
        # Triangles.
        'tri1': {
            'tg': stk.cage.M3L3Triangle(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    3: 0, 4: 0, 5: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        'tri2': {
            'tg': stk.cage.M3L3Triangle(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    3: 0, 4: 0, 5: 1,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        # Squares.
        'sqr1': {
            'tg': stk.cage.M4L4Square(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    4: 0, 5: 0, 6: 0, 7: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr2': {
            'tg': stk.cage.M4L4Square(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    4: 0, 5: 1, 6: 0, 7: 1,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr3': {
            'tg': stk.cage.M4L4Square(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    4: 0, 5: 0, 6: 1, 7: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr4': {
            'tg': stk.cage.M4L4Square(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    4: 0, 5: 0, 6: 1, 7: 1,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        # Hexagons.
        'hex1': {
            'tg': M6L6(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    6: 0, 7: 0, 8: 0, 9: 0, 10: 0, 11: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex2': {
            'tg': M6L6(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    6: 0, 7: 1, 8: 0, 9: 0, 10: 0, 11: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex3': {
            'tg': M6L6(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    6: 0, 7: 1, 8: 1, 9: 0, 10: 0, 11: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex4': {
            'tg': M6L6(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    6: 0, 7: 1, 8: 0, 9: 1, 10: 0, 11: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex5': {
            'tg': M6L6(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    6: 0, 7: 1, 8: 1, 9: 1, 10: 0, 11: 0,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex6': {
            'tg': M6L6(
                corners=corner_bb,
                linkers=lig_bb,
                vertex_alignments={
                    6: 0, 7: 1, 8: 0, 9: 1, 10: 0, 11: 1,
                },
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
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
