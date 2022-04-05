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
import networkx as nx
import pymatgen.core as pmg
from pymatgen.analysis.local_env import (
    LocalStructOrderParams,
)
import json

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


def convert_stk_to_pymatgen(stk_mol):
    """
    Convert stk.Molecule to pymatgen.Molecule.

    Parameters
    ----------
    stk_mol : :class:`stk.Molecule`
        Stk molecule to convert.

    Returns
    -------
    pmg_mol : :class:`pymatgen.Molecule`
        Corresponding pymatgen Molecule.

    """
    stk_mol.write('temp.xyz')
    pmg_mol = pmg.Molecule.from_file('temp.xyz')
    os.system('rm temp.xyz')

    return pmg_mol


def calculate_sites_order_values(
    molecule,
    site_idxs,
    target_species_type=None,
    neigh_idxs=None
):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    molecule : :class:`pmg.Molecule` or :class:`pmg.Structure`
        Pymatgen (pmg) molecule/structure to analyse.

    site_idxs : :class:`list` of :class:`int`
        Atom ids of sites to calculate OP of.

    target_species_type : :class:`str`
        Target neighbour element to use in OP calculation.
        Defaults to :class:`NoneType` if no target species is known.

    neigh_idxs : :class:`list` of :class:`list` of :class:`int`
        Neighbours of each atom in site_idx. Ordering is important.
        Defaults to :class:`NoneType` for when using
        :class:`pmg.Structure` - i.e. a structure with a lattice.

    Returns
    -------
    results : :class:`dict`
        Dictionary of format
        site_idx: dict of order parameters
        {
            `oct`: :class:`float`,
            `sq_plan`: :class:`float`,
            `q2`: :class:`float`,
            `q4`: :class:`float`,
            `q6`: :class:`float`
        }.

    """

    results = {}

    if target_species_type is None:
        targ_species = None
    else:
        targ_species = pmg.Species(target_species_type)

    # Define local order parameters class based on desired types.
    types = [
        'sq_plan',  # Square planar envs.
    ]
    loc_ops = LocalStructOrderParams(
        types=types,
    )
    if neigh_idxs is None:
        for site in site_idxs:
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                target_spec=[targ_species]
            )
            results[site] = {i: j for i, j in zip(types, site_results)}
    else:
        for site, neigh in zip(site_idxs, neigh_idxs):
            site_results = loc_ops.get_order_parameters(
                structure=molecule,
                n=site,
                indices_neighs=neigh,
                target_spec=targ_species
            )
            results[site] = {i: j for i, j in zip(types, site_results)}

    return results


def get_order_values(mol, metal, per_site=False):
    """
    Calculate order parameters around metal centres.

    Parameters
    ----------
    mol : :class:`stk.ConstructedMolecule`
        stk molecule to analyse.

    metal : :class:`int`
        Element number of metal atom.

    per_site : :class:`bool`
        Defaults to False. True if the OPs for each site are desired.

    Returns
    -------
    results : :class:`dict`
        Dictionary of order parameter max/mins/averages if `per_site`
        is False.

    """

    pmg_mol = convert_stk_to_pymatgen(stk_mol=mol)
    # Get sites of interest and their neighbours.
    sites = []
    neighs = []
    for atom in mol.get_atoms():
        if atom.get_atomic_number() == metal:
            sites.append(atom.get_id())
            bonds = [
                i
                for i in mol.get_bonds()
                if i.get_atom1().get_id() == atom.get_id()
                or i.get_atom2().get_id() == atom.get_id()
            ]
            a_neigh = []
            for b in bonds:
                if b.get_atom1().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom2().get_id())
                elif b.get_atom2().get_id() == atom.get_id():
                    a_neigh.append(b.get_atom1().get_id())
            neighs.append(a_neigh)

    order_values = calculate_sites_order_values(
        molecule=pmg_mol,
        site_idxs=sites,
        neigh_idxs=neighs,
    )

    if per_site:
        results = order_values
        return results
    else:
        # Get max, mins and averages of all OPs for the whole molecule.
        OPs = [order_values[i].keys() for i in order_values][0]
        OP_lists = {}
        for OP in OPs:
            OP_lists[OP] = [order_values[i][OP] for i in order_values]

        results = {
            # OP: (min, max, avg)
            i: {
                'min': min(OP_lists[i]),
                'max': max(OP_lists[i]),
                'avg': np.average(OP_lists[i])
            }
            for i in OP_lists
        }

        return results


def get_organic_linkers(
    cage,
    metal_atom_nos,
    calc_dir,
    file_prefix=None,
):
    """
    Extract a list of organic linker .Molecules from a cage.

    Parameters
    ----------
    cage : :class:`stk.Molecule`
        Molecule to get the organic linkers from.

    metal_atom_nos : :class:`iterable` of :class:`int`
        The atomic number of metal atoms to remove from structure.

    file_prefix : :class:`str`, optional
        Prefix to file name of each output ligand structure.
        Eventual file name is:
        "file_prefix"{number of atoms}_{idx}_{i}.mol
        Where `idx` determines if a molecule is unique by smiles.

    Returns
    -------
    org_lig : :class:`dict` of :class:`stk.BuildingBlock`
        Dictionary of building blocks where the key is the file name,
        and the value is the stk building block.

    smiles_keys : :class:`dict` of :class:`int`
        Key is the linker smiles, value is the idx of that smiles.

    """

    org_lig = {}

    # Produce a graph from the cage that does not include metals.
    cage_g = nx.Graph()
    atom_ids_in_G = set()
    for atom in cage.get_atoms():
        if atom.get_atomic_number() in metal_atom_nos:
            continue
        cage_g.add_node(atom)
        atom_ids_in_G.add(atom.get_id())

    # Add edges.
    for bond in cage.get_bonds():
        a1id = bond.get_atom1().get_id()
        a2id = bond.get_atom2().get_id()
        if a1id in atom_ids_in_G and a2id in atom_ids_in_G:
            cage_g.add_edge(bond.get_atom1(), bond.get_atom2())

    # Get disconnected subgraphs as molecules.
    # Sort and sort atom ids to ensure molecules are read by RDKIT
    # correctly.
    connected_graphs = [
        sorted(subgraph, key=lambda a: a.get_id())
        for subgraph in sorted(nx.connected_components(cage_g))
    ]
    smiles_keys = {}
    for i, cg in enumerate(connected_graphs):
        # Get atoms from nodes.
        atoms = list(cg)
        atom_ids = [i.get_id() for i in atoms]
        cage.write(
            'temporary_linker.mol',
            atom_ids=atom_ids,
        )
        temporary_linker = stk.BuildingBlock.init_from_file(
            'temporary_linker.mol'
        ).with_canonical_atom_ordering()
        smiles_key = stk.Smiles().get_key(temporary_linker)
        if smiles_key not in smiles_keys:
            smiles_keys[smiles_key] = len(smiles_keys.values())+1
        idx = smiles_keys[smiles_key]
        sgt = str(len(atoms))
        # Write to mol file.
        if file_prefix is None:
            filename_ = f'organic_linker_s{sgt}_{idx}_{i}.mol'
        else:
            filename_ = f'{file_prefix}{sgt}_{idx}_{i}.mol'

        org_lig[filename_] = temporary_linker
        os.system('rm temporary_linker.mol')
        # Rewrite to fix atom ids.
        org_lig[filename_].write(os.path.join(calc_dir, filename_))
        org_lig[filename_] = stk.BuildingBlock.init_from_file(
            path=os.path.join(calc_dir, filename_)
        )

    return org_lig, smiles_keys


def calculate_ligand_SE(
    org_ligs,
    lowe_ligand_energy,
    output_json,
    calc_dir,
):

    # Check if output file exists.
    if not os.path.exists(output_json):
        logging.info(
            f'calculating strain energy of target with 30 atoms.'
        )
        strain_energies = {}
        # Iterate over ligands.
        for lig in org_ligs:
            stk_lig = org_ligs[lig]
            # Only run for the target ligand.
            if stk_lig.get_num_atoms() != 30:
                continue

            # Calculate energy of extracted ligand.
            energy_au = get_energy(
                molecule=stk_lig,
                name=lig,
                charge=0,
                calc_dir=calc_dir,
            )
            # kJ/mol.
            E_extracted = energy_au * 2625.5
            E_free = lowe_ligand_energy * 2625.5
            # Add to list the strain energy:
            # (E(extracted) - E(optimised/free))
            lse = E_extracted - E_free
            strain_energies[lig] = lse

        # Write data.
        with open(output_json, 'w') as f:
            json.dump(strain_energies, f)

    # Get data.
    with open(output_json, 'r') as f:
        strain_energies = json.load(f)

    return strain_energies


def get_energy(molecule, name, charge, calc_dir):
    output_dir = os.path.join(calc_dir, f'{name}_xtbey')
    output_file = os.path.join(calc_dir, f'{name}_xtb.ey')
    if os.path.exists(output_file):
        with open(output_file, 'r') as f:
            lines = f.readlines()
        for line in lines:
            energy = float(line.rstrip())
            break
    else:
        logging.info(f'xtb energy calculation of {name}')
        xtb = stko.XTBEnergy(
            xtb_path=xtb_path(),
            output_dir=output_dir,
            gfn_version=2,
            num_cores=6,
            charge=charge,
            num_unpaired_electrons=0,
            unlimited_memory=True,
        )
        energy = xtb.get_energy(mol=molecule)
        with open(output_file, 'w') as f:
            f.write(f'{energy}\n')

    # In a.u.
    return energy
