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
from utilities import AromaticCNCFactory, get_lowest_energy_conformer


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


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
