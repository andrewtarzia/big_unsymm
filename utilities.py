#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for utility functions.

Author: Andrew Tarzia

"""

import logging
import os
import stk
import stko
import numpy as np

from env_set import xtb_path, crest_path


class AromaticCNCFactory(stk.FunctionalGroupFactory):
    """
    A subclass of stk.SmartsFunctionalGroupFactory.

    """

    def __init__(self, bonders=(1, ), deleters=()):
        """
        Initialise :class:`.AromaticCNCFactory`.

        """

        self._bonders = bonders
        self._deleters = deleters

    def get_functional_groups(self, molecule):
        generic_functional_groups = stk.SmartsFunctionalGroupFactory(
            smarts='[#6]~[#7X2]~[#6]',
            bonders=self._bonders,
            deleters=self._deleters,
        ).get_functional_groups(molecule)
        for fg in generic_functional_groups:
            atom_ids = (i.get_id() for i in fg.get_atoms())
            atoms = tuple(molecule.get_atoms(atom_ids))
            yield AromaticCNC(
                carbon1=atoms[0],
                nitrogen=atoms[1],
                carbon2=atoms[2],
                bonders=tuple(atoms[i] for i in self._bonders),
                deleters=tuple(atoms[i] for i in self._deleters),
            )


class AromaticCNC(stk.GenericFunctionalGroup):
    """
    Represents an N atom in pyridine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon][nitrogen][carbon]``.

    """

    def __init__(
        self,
        carbon1,
        nitrogen,
        carbon2,
        bonders,
        deleters,
    ):
        """
        Initialize a :class:`.AromaticCNC` instance.

        Parameters
        ----------
        carbon1 : :class:`.C`
            The first carbon atom.

        nitrogen : :class:`.N`
            The nitrogen atom.

        carbon2 : :class:`.C`
            The second carbon atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon1 = carbon1
        self._nitrogen = nitrogen
        self._carbon2 = carbon2
        atoms = (carbon1, nitrogen, carbon2)
        super().__init__(atoms, bonders, deleters)

    def get_carbon1(self):
        return self._carbon1

    def get_carbon2(self):
        return self._carbon2

    def get_nitrogen(self):
        return self._nitrogen

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._nitrogen = self._nitrogen
        clone._carbon2 = self._carbon2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._nitrogen}, {self._carbon2}, '
            f'bonders={self._bonders})'
        )


class MissingSettingError(Exception):
    ...


def get_lowest_energy_conformer(name, mol, settings, calc_dir):
    """
    Get lowest energy conformer of molecule.

    Method:
        1) squick CREST conformer search
        2) xTB `opt_level` optimisation of lowest energy conformer
        3) save file

    """

    # Check for missing settings.
    req_settings = [
        'final_opt_level', 'conf_opt_level', 'charge', 'no_unpaired_e',
        'max_runs', 'calc_hessian', 'solvent', 'nc',
        'etemp', 'keepdir', 'cross', 'md_len', 'ewin', 'speed_setting'
    ]
    for i in req_settings:
        if i not in settings:
            raise MissingSettingError(
                f'Settings missing {i}. Has {settings.keys()}.'
            )

    crest_output_dir = os.path.join(calc_dir, f'{name}_crest')
    crest_opt_dir = os.path.join(calc_dir, f'{name}_crestxtb')
    low_unopt_file = os.path.join(calc_dir, f'{name}_low_e_unopt.mol')
    low_opt_file = os.path.join(calc_dir, f'{name}_low_e_opt.mol')

    if not os.path.exists(low_unopt_file):
        logging.info(f'running crest opt of {name}')
        xtb_crest = stko.XTBCREST(
            crest_path=crest_path(),
            xtb_path=xtb_path(),
            gfn_version=2,
            output_dir=crest_output_dir,
            num_cores=settings['nc'],
            ewin=settings['ewin'],
            opt_level=settings['conf_opt_level'],
            charge=settings['charge'],
            electronic_temperature=settings['etemp'],
            num_unpaired_electrons=settings['no_unpaired_e'],
            solvent=settings['solvent'],
            keepdir=settings['keepdir'],
            cross=settings['cross'],
            speed_setting=settings['speed_setting'],
            md_len=settings['md_len'],
            unlimited_memory=True,
        )
        low_e_conf = xtb_crest.optimize(mol=mol)
        # Save lowest energy conformer.
        low_e_conf.write(low_unopt_file)
    else:
        low_e_conf = stk.BuildingBlock.init_from_file(low_unopt_file)

    if not os.path.exists(low_opt_file):
        logging.info(f'running xtb opt of crest conformer of {name}')
        xtb_opt = stko.XTB(
            xtb_path=xtb_path(),
            output_dir=crest_opt_dir,
            gfn_version=2,
            num_cores=6,
            opt_level=settings['final_opt_level'],
            charge=settings['charge'],
            num_unpaired_electrons=settings['no_unpaired_e'],
            max_runs=settings['max_runs'],
            calculate_hessian=settings['calc_hessian'],
            unlimited_memory=True,
            solvent=settings['solvent'],
        )
        low_e_opt = xtb_opt.optimize(mol=low_e_conf)
        low_e_opt.write(low_opt_file)
    else:
        low_e_opt = stk.BuildingBlock.init_from_file(low_opt_file)

    return low_e_opt


def split_xyz_file(num_atoms, xyz_file):
    """
    Splits xyz trajectory file into xyz files.

    """

    with open(xyz_file, 'r') as f:
        lines = f.readlines()

    file_strings = []
    string = []
    for line in lines:
        if f' {num_atoms} ' in f' {line.strip()} ':
            if len(string) == 0:
                string.append(line)
            else:
                # New block.
                file_strings.append(string)
                string = [line]
        else:
            string.append(line)
    # Add last set.
    file_strings.append(string)

    out_files = []
    for i, fs in enumerate(file_strings):
        file_name = xyz_file.replace('.xyz', f'_s{i}.xyz')
        with open(file_name, 'w') as f:
            for line in fs:
                f.write(line)
        out_files.append(file_name)

    return out_files


def unit_vector(vector):
    """
    Returns the unit vector of the vector.

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    """
    return vector / np.linalg.norm(vector)


def angle_between(v1, v2, normal=None):
    """
    Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

    https://stackoverflow.com/questions/2827393/
    angles-between-two-n-dimensional-vectors-in-python/
    13849249#13849249

    If normal is given, the angle polarity is determined using the
    cross product of the two vectors.

    """

    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    if normal is not None:
        # Get normal vector and cross product to determine sign.
        cross = np.cross(v1_u, v2_u)
        if np.dot(normal, cross) < 0:
            angle = -angle
    return angle
