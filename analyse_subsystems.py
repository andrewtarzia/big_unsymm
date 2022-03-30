#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse all subsystems constructed.

Author: Andrew Tarzia

"""

import logging
import glob
import sys
import os
import json
import stk

from env_set import calc_path, liga_path, subs_path
from utilities import (
    get_order_values,
    get_energy,
)
from plotting import plot_energies, plot_ops


def get_min_order_parameter(molecule):

    order_results = get_order_values(
        mol=molecule,
        metal=46
    )
    return order_results['sq_plan']['min']


def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 1 arguments:'
        )
        sys.exit()
    else:
        pass

    li_path = liga_path()
    lowe_ligand = stk.BuildingBlock.init_from_file(
        path=os.path.join(li_path, 'lig_lowe.mol'),
    )
    _wd = subs_path()
    _cd = calc_path()

    structure_files = glob.glob(os.path.join(_wd, '*_opt.mol'))
    logging.info(f'there are {len(structure_files)} structures.')
    structure_results = {
        i.split('/')[-1].replace('_opt.mol', ''): {}
        for i in structure_files
    }
    structure_res_file = os.path.join(_wd, 'all_structure_res.json')
    if os.path.exists(structure_res_file):
        with open(structure_res_file, 'r') as f:
            structure_results = json.load(f)
    else:
        for s_file in structure_files:
            name = s_file.split('/')[-1].replace('_opt.mol', '')
            molecule = stk.BuildingBlock.init_from_file(s_file)

            if 'tri' in name:
                charge = 6
                exp_lig = None
            elif 'sqr' in name:
                charge = 8
                exp_lig = None
            elif 'hex' in name:
                charge = 12
                exp_lig = None
                continue

            xtb_energy = get_energy(molecule, name, charge, _cd)
            structure_results[name]['xtb_energy'] = xtb_energy

            min_order_param = get_min_order_parameter(molecule)
            structure_results[name]['min_order_param'] = (
                min_order_param
            )

        with open(structure_res_file, 'w') as f:
            json.dump(structure_results, f)

    print(structure_results)

    plot_energies(
        results_dict=structure_results,
        outname='subs_energies',
        per_ligand=True,
    )
    plot_ops(
        results_dict=structure_results,
        outname='subs_ops',
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
