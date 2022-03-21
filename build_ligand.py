#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build the ligand in this project.

Author: Andrew Tarzia

"""

import logging
import sys
import os
import stk
import numpy as np

from env_set import liga_path, calc_path, crest_conformer_settings
from utilities import (
    angle_between,
    AromaticCNCFactory,
    AromaticCNC,
    get_lowest_energy_conformer,
    split_xyz_file,
)


def calculate_N_COM_N_angle(bb):
    """
    Calculate the N-COM-N angle of a ditopic building block.

    This function will not work for cages built from FGs other than
    metals + AromaticCNC.

    Parameters
    ----------
    bb : :class:`stk.BuildingBlock`
        stk molecule to analyse.

    Returns
    -------
    angle : :class:`float`
        Angle between two bonding vectors of molecule.

    """

    fg_counts = 0
    fg_positions = []
    for fg in bb.get_functional_groups():
        if isinstance(fg, AromaticCNC):
            fg_counts += 1
            # Get geometrical properties of the FG.
            # Get N position - deleter.
            N_position, = bb.get_atomic_positions(
                atom_ids=fg.get_nitrogen().get_id()
            )
            fg_positions.append(N_position)

    if fg_counts != 2:
        raise ValueError(
            f'{bb} does not have 2 AromaticCNC or AromaticCNN '
            'functional groups.'
        )

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
    logging.info(
        f'getting optimal conformer of {name} from {crest_output_dir}'
    )

    # Analyse all conformers from CREST.
    crest_conformer_files = split_xyz_file(
        num_atoms=mol.get_num_atoms(),
        xyz_file=(
            os.path.join(crest_output_dir, 'crest_conformers.xyz'),
        ),
    )
    logging.info(f'{name} has {len(crest_conformer_files)} conformers')

    min_angle = calculate_N_COM_N_angle(mol)
    for cre_file in crest_conformer_files:
        logging.info(cre_file)
        _temp_mol = mol.with_structure_from_file(cre_file)
        angle = calculate_N_COM_N_angle(_temp_mol)
        if angle < min_angle:
            min_angle = angle
            out_molecule = mol.with_structure_from_file(cre_file)

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

    _wd = liga_path()
    _cd = calc_path()
    _crest_settings = crest_conformer_settings()

    if not os.path.exists(_wd):
        os.mkdir(_wd)

    if not os.path.exists(_cd):
        os.mkdir(_cd)

    ligand_smiles = 'C1=CC(=CN=C1)C2=CC=C(C=C2)C3=CC=NC=C3'
    unopt_mol = stk.BuildingBlock(
        smiles=ligand_smiles,
        functional_groups=(AromaticCNCFactory(), ),
    )
    unopt_mol.write(os.path.join(_wd, f'lig_unopt.mol'))

    logging.info(f'running crest conformer search: {_crest_settings}')
    lowe_mol = get_lowest_energy_conformer(
        name='lig',
        mol=unopt_mol,
        settings=_crest_settings,
        calc_dir=_cd,
    )
    lowe_mol.write(os.path.join(_wd, f'lig_lowe.mol'))

    chosen_mol = get_chosen_conformer(
        mol=unopt_mol,
        name='lig',
        calc_dir=_cd,
    )
    chosen_mol.write(os.path.join(_wd, f'lig_chosen.mol'))


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
