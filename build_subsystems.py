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
import numpy as np

from env_set import subs_path, calc_path, meta_path
from optimisation import subsystem_optimisation_sequence
from utilities import AromaticCNC


def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 0 arguments:'
        )
        sys.exit()
    else:
        pass

    me_path = meta_path()
    pm_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 'm1_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )
    mp_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 'm4_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )
    mm_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 'm2_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )
    pp_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 'm3_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )

    tpm_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 't1_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )
    tmp_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 't4_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )
    tmm_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 't2_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )
    tpp_bb = stk.BuildingBlock.init_from_file(
        path=os.path.join(me_path, 't3_opt.mol'),
        functional_groups=(stk.BromoFactory(), ),
    )

    _wd = subs_path()
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
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(pm_bb, mm_bb, pp_bb),
                repeating_unit='AAA',
                num_repeating_units=1,
                orientations=(0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        'tri2': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(pm_bb, mm_bb, pp_bb),
                repeating_unit='ABC',
                num_repeating_units=1,
                orientations=(0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*3,
        },
        # Squares.
        'sqr1': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(pm_bb, mm_bb, pp_bb, mp_bb),
                repeating_unit='AAAA',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr2': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(pm_bb, mm_bb, pp_bb, mp_bb),
                repeating_unit='BCBC',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr3': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(pm_bb, mm_bb, pp_bb, mp_bb),
                repeating_unit='AABC',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        'sqr4': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(pm_bb, mm_bb, pp_bb, mp_bb),
                repeating_unit='BDCA',
                num_repeating_units=1,
                orientations=(0, 1, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*4,
        },
        # Hexagons.
        'hex1': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(tpm_bb, tmm_bb, tpp_bb, tmp_bb),
                repeating_unit='AAAAAA',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex2': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(tpm_bb, tmm_bb, tpp_bb, tmp_bb),
                repeating_unit='ACBAAA',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex3': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(tpm_bb, tmm_bb, tpp_bb, tmp_bb),
                repeating_unit='ACDBAA',
                num_repeating_units=1,
                orientations=(0, 0, 1, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex4': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(tpm_bb, tmm_bb, tpp_bb, tmp_bb),
                repeating_unit='ACBCBA',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex5': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(tpm_bb, tmm_bb, tpp_bb, tmp_bb),
                repeating_unit='ACDDBA',
                num_repeating_units=1,
                orientations=(0, 0, 1, 1, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
        },
        'hex6': {
            'tg': stk.macrocycle.Macrocycle(
                building_blocks=(tpm_bb, tmm_bb, tpp_bb, tmp_bb),
                repeating_unit='BCBCBC',
                num_repeating_units=1,
                orientations=(0, 0, 0, 0, 0, 0),
                reaction_factory=_react_factory,
                optimizer=_init_opts,
            ),
            'charge': 2*6,
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

        if not os.path.exists(opt_file):
            logging.info(f'optimising {topo}')
            opt_mol = subsystem_optimisation_sequence(
                mol=unopt_mol,
                name=topo,
                charge=charge,
                calc_dir=_cd,
            )
            opt_mol = opt_mol.with_centroid(np.array((0, 0, 0)))
            opt_mol.write(opt_file)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
