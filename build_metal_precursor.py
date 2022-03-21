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
import numpy as np

from env_set import (
    meta_path,
    calc_path,
    xtb_path,
    crest_conformer_settings,
)
from utilities import (
    AromaticCNCFactory,
    angle_between,
    split_xyz_file,
    get_lowest_energy_conformer,
)


def calculate_Br_COM_Br_angle(bb):
    fg_counts = 0
    fg_positions = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, stk.Bromo):
            fg_counts += 1
            Br_position, = bb.get_atomic_positions(
                atom_ids=fg.get_bromine().get_id()
            )
            fg_positions.append(Br_position)

    if fg_counts != 2:
        raise ValueError(f'{bb} does not have 2 bromines.')

    # Get building block centroid.
    centroid = bb.get_centroid()

    # Get vectors.
    fg_vectors = [i-centroid for i in fg_positions]

    # Calculate the angle between the two vectors.
    angle = np.degrees(angle_between(*fg_vectors))
    return angle


def get_chosen_conformer(mol, name, calc_dir):
    """
    Select and optimize a conformer with desired directionality.

    Currently:
        Best directionality will be defined by the smallest
        N-ligand centroid-N angle.

    """

    crest_output_dir = os.path.join(calc_dir, f'{name}_crest')
    _crest_settings = crest_conformer_settings()

    settings = _crest_settings
    settings['speed_setting'] = 'mquick'
    settings['charge'] = 2
    settings['final_opt_level'] = 'normal'
    logging.info(f'running crest conformer search: {settings}')
    lowe_mol = get_lowest_energy_conformer(
        name=name,
        mol=mol,
        settings=settings,
        calc_dir=calc_dir,
    )

    logging.info(
        f'getting optimal conformer of {name} from {crest_output_dir}'
    )
    lowe_mol = stk.BuildingBlock.init_from_molecule(
        molecule=lowe_mol,
        functional_groups=(stk.BromoFactory(), ),
    )
    # Analyse all conformers from CREST.
    crest_conformer_files = split_xyz_file(
        num_atoms=lowe_mol.get_num_atoms(),
        xyz_file=(
            os.path.join(crest_output_dir, 'crest_conformers.xyz')
        ),
    )
    logging.info(f'{name} has {len(crest_conformer_files)} conformers')

    min_angle = calculate_Br_COM_Br_angle(lowe_mol)
    print(min_angle)
    for cre_file in crest_conformer_files:
        _temp_mol = lowe_mol.with_structure_from_file(cre_file)
        angle = calculate_Br_COM_Br_angle(_temp_mol)
        print(angle)
        if angle < min_angle:
            min_angle = angle
            out_molecule = lowe_mol.with_structure_from_file(cre_file)

    return out_molecule


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

    pyridine = stk.BuildingBlock(
        smiles='C1=CC=NC=C1',
        functional_groups=(AromaticCNCFactory(), ),
    )
    meta = stk.BuildingBlock(
        smiles='C1=CC(=CN=C1)Br',
        functional_groups=(AromaticCNCFactory(), ),
    )
    para = stk.BuildingBlock(
        smiles='C1=CN=CC=C1Br',
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

    _topos = {
        'cis1m1p': {
            pyridine: (0, 1),
            para: (2, ),
            meta: (3, ),
        },
        'cis2m0p': {
            pyridine: (0, 1),
            meta: (2, 3),
        },
        'cis0m2p': {
            pyridine: (0, 1),
            para: (2, 3),
        },
        'trans1m1p': {
            pyridine: (0, 2),
            para: (1, ),
            meta: (3, ),
        },
        'trans2m0p': {
            pyridine: (0, 2),
            meta: (1, 3),
        },
        'trans0m2p': {
            pyridine: (0, 2),
            para: (1, 3),
        },
    }

    for topo in _topos:
        meta_unopt = os.path.join(_wd, f'{topo}_unopt.mol')
        meta_opt = os.path.join(_wd, f'{topo}_opt.mol')
        meta_fin = os.path.join(_wd, f'{topo}_final.mol')

        sqpl_unopt = stk.ConstructedMolecule(
            topology_graph=stk.metal_complex.SquarePlanar(
                metals=pd,
                ligands=_topos[topo],
                optimizer=stk.MCHammer(
                    num_steps=300,
                    target_bond_length=1.8,
                ),
            ),
        )
        sqpl_unopt.write(meta_unopt)

        if not os.path.exists(meta_opt):
            logging.info(f'xtb optimisation of {topo}')
            xtb_opt = stko.XTB(
                xtb_path=xtb_path(),
                output_dir=os.path.join(_cd, f'{topo}_xtbopt'),
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
        else:
            sqpl_opt = sqpl_unopt.with_structure_from_file(meta_opt)

        if '2m' in topo:
            if os.path.exists(meta_fin):
                continue
            sqpl_fin = get_chosen_conformer(
                mol=sqpl_opt,
                name=topo,
                calc_dir=_cd,
            )
            sqpl_fin = xtb_opt.optimize(mol=sqpl_opt)
            sqpl_fin.write(meta_fin)
        else:
            sqpl_opt.write(meta_fin)


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
