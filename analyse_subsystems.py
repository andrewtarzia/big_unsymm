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
    get_organic_linkers,
    calculate_ligand_SE,
)
from plotting import plot_subs_energies, plot_ops, plot_strain_energies


def get_min_order_parameter(molecule):

    order_results = get_order_values(
        mol=molecule,
        metal=46
    )
    return order_results['sq_plan']['min']



class UnexpectedNumLigands(Exception):
    ...


def get_subs_sum_strain_energy(
    molecule,
    name,
    exp_lig,
    lowe_ligand,
    calc_dir,
):

    ls_file = os.path.join(calc_dir, f'{name}_strain.json')

    org_ligs, smiles_keys = get_organic_linkers(
        cage=molecule,
        metal_atom_nos=(46, ),
        file_prefix=f'{name}_sg',
        calc_dir=calc_dir,
    )

    num_unique_ligands = len(set(smiles_keys.values()))
    if num_unique_ligands != exp_lig:
        raise UnexpectedNumLigands(
            f'{name} had {num_unique_ligands} unique ligands'
            f', {exp_lig} were expected. Suggests bad '
            'optimization. Recommend reoptimising structure.'
        )


    lowe_ligand_energy = get_energy(
        molecule=lowe_ligand,
        name='lig',
        charge=0,
        calc_dir=calc_dir,
    )

    lse_dict = calculate_ligand_SE(
        org_ligs=org_ligs,
        lowe_ligand_energy=lowe_ligand_energy,
        output_json=f'{ls_file}.json',
        calc_dir=calc_dir,
    )

    sum_strain = sum(lse_dict.values())
    return sum_strain


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

    structure_files = sorted(glob.glob(os.path.join(_wd, '*_opt.mol')))
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
            elif 'sqr' in name:
                charge = 8
            elif 'hex' in name:
                charge = 12

            xtb_energy = get_energy(molecule, name, charge, _cd)
            structure_results[name]['xtb_energy'] = xtb_energy

            min_order_param = get_min_order_parameter(molecule)
            structure_results[name]['min_order_param'] = (
                min_order_param
            )

            sum_strain_energy = get_subs_sum_strain_energy(
                molecule=molecule,
                name=name,
                exp_lig=2,
                lowe_ligand=lowe_ligand,
                calc_dir=_cd,
            )
            structure_results[name]['sum_strain_energy'] = (
                sum_strain_energy
            )

        with open(structure_res_file, 'w') as f:
            json.dump(structure_results, f)

    print(structure_results)

    plot_subs_energies(
        results_dict=structure_results,
        outname='subs_energies',
    )
    plot_ops(
        results_dict=structure_results,
        outname='subs_ops',
    )
    plot_strain_energies(
        results_dict=structure_results,
        outname='subs_strain_energies',
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
