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

from env_set import cage_path, calc_path, meta_path
from optimisation import optimisation_sequence
from m6l6 import M6L6


def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 0 arguments:'
        )
        sys.exit()
    else:
        pass

    lig_bb = stk.BuildingBlock(
        smiles='C1=CC(=CC=C1Br)Br',
        functional_groups=(stk.BromoFactory(), ),
    )
    me_path = meta_path()
    corner_bbs = {
        'cis1m1p': stk.BuildingBlock.init_from_file(
            path=os.path.join(me_path, 'cis1m1p_opt.mol'),
            functional_groups=(stk.BromoFactory(), ),
        ),
        'cis2m0p': stk.BuildingBlock.init_from_file(
            path=os.path.join(me_path, 'cis2m0p_opt.mol'),
            functional_groups=(stk.BromoFactory(), ),
        ),
        'cis0m2p': stk.BuildingBlock.init_from_file(
            path=os.path.join(me_path, 'cis0m2p_opt.mol'),
            functional_groups=(stk.BromoFactory(), ),
        ),
        'trans1m1p': stk.BuildingBlock.init_from_file(
            path=os.path.join(me_path, 'trans1m1p_opt.mol'),
            functional_groups=(stk.BromoFactory(), ),
        ),
        'trans2m0p': stk.BuildingBlock.init_from_file(
            path=os.path.join(me_path, 'trans2m0p_opt.mol'),
            functional_groups=(stk.BromoFactory(), ),
        ),
        'trans0m2p': stk.BuildingBlock.init_from_file(
            path=os.path.join(me_path, 'trans0m2p_opt.mol'),
            functional_groups=(stk.BromoFactory(), ),
        ),
    }
    _wd = cage_path()
    _cd = calc_path()

    if not os.path.exists(_wd):
        os.mkdir(_wd)

    if not os.path.exists(_cd):
        os.mkdir(_cd)

    # Define series of topologies to build.
    _init_opts = stk.NullOptimizer()
    # _init_opts = stk.MCHammer(target_bond_length=2.5, num_steps=1000)
    _topos = {
        # Triangles.
        'tri1': {
            'tg': stk.cage.M3L3Triangle(
                corners={
                    corner_bbs['cis0m2p']: (0, 1, 2),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        'tri2': {
            'tg': stk.cage.M3L3Triangle(
                corners={
                    corner_bbs['cis2m0p']: (0, 1, 2),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        'tri3': {
            'tg': stk.cage.M3L3Triangle(
                corners={
                    corner_bbs['cis1m1p']: (0, 1, 2),
                },
                linkers=lig_bb,
                vertex_alignments={
                    0: 1,
                },
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        'tri4': {
            'tg': stk.cage.M3L3Triangle(
                corners={
                    corner_bbs['cis0m2p']: (0, ),
                    corner_bbs['cis2m0p']: (1, ),
                    corner_bbs['cis1m1p']: (2, ),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        # Squares.
        'sqr1': {
            'tg': stk.cage.M4L4Square(
                corners={
                    corner_bbs['cis0m2p']: (0, 1, 2, 3),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr2': {
            'tg': stk.cage.M4L4Square(
                corners={
                    corner_bbs['cis2m0p']: (0, 1, 2, 3),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr3': {
            'tg': stk.cage.M4L4Square(
                corners={
                    corner_bbs['cis2m0p']: (0, 2),
                    corner_bbs['cis0m2p']: (1, 3),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr4': {
            'tg': stk.cage.M4L4Square(
                corners={
                    corner_bbs['cis1m1p']: (0, 1, 2, 3),
                },
                linkers=lig_bb,
                vertex_alignments={0: 1},
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr5': {
            'tg': stk.cage.M4L4Square(
                corners={
                    corner_bbs['cis2m0p']: (0, ),
                    corner_bbs['cis0m2p']: (1, ),
                    corner_bbs['cis1m1p']: (2, 3),
                },
                linkers=lig_bb,
                vertex_alignments={2: 1, 3: 1},
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr6': {
            'tg': stk.cage.M4L4Square(
                corners={
                    corner_bbs['cis2m0p']: (0, ),
                    corner_bbs['cis0m2p']: (2, ),
                    corner_bbs['cis1m1p']: (1, 3),
                },
                linkers=lig_bb,
                vertex_alignments={3: 1},
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        # Hexagons.
        'hex1': {
            'tg': M6L6(
                corners={
                    corner_bbs['trans0m2p']: (0, 1, 2, 3, 4, 5),
                },
                linkers=lig_bb,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex2': {
            'tg': M6L6(
                corners={
                    corner_bbs['trans2m0p']: (0, 1, 2, 3, 4, 5),
                },
                linkers=lig_bb,
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
