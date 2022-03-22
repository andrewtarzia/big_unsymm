#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for spindry functions.

Author: Andrew Tarzia

"""

import spindry as spd
import numpy as np

from utilities import angle_between


class FGPotential(spd.Potential):

    def __init__(self, moving_pairs, target):
        self._moving_pairs = moving_pairs
        self._target = target

    def compute_potential(self, supramolecule):
        pos_mat = supramolecule.get_position_matrix()

        potential = 0
        for pair in self._moving_pairs:
            p1_pos = pos_mat[pair[1]]
            p2_pos = pos_mat[pair[0]]
            pair_v = p2_pos - p1_pos
            # Calculate the angle between the two vectors.
            angle = np.degrees(angle_between(self._target, pair_v))
            potential += angle

        return potential


def get_supramolecule(mol):

    # Remove Pd - N bonds.
    bonds = tuple(
        i for i in mol.get_bonds()
        if i.get_order() != 9
    )

    supramolecule = spd.SupraMolecule(
        atoms=(
            spd.Atom(
                id=atom.get_id(),
                element_string=atom.__class__.__name__,
            ) for atom in mol.get_atoms()
        ),
        bonds=(
            spd.Bond(
                id=i,
                atom_ids=(
                    bond.get_atom1().get_id(),
                    bond.get_atom2().get_id(),
                )
            ) for i, bond in enumerate(bonds)
        ),
        position_matrix=mol.get_position_matrix(),
    )
    return supramolecule


def rotate_fgs(mol):
    smolecule = get_supramolecule(mol)

    # Get rotatable components.
    bromo_ids = tuple(
        fg.get_bromine().get_id()
        for fg in mol.get_functional_groups()
    )
    movable_components = {}
    for i, comp in enumerate(smolecule.get_components()):
        comp_atom_ids = set(i.get_id() for i in comp.get_atoms())
        movable = False
        for b_id in bromo_ids:
            if b_id in comp_atom_ids:
                movable = True
        movable_components[i] = movable

    movable_components = tuple(
        i for i in movable_components if movable_components[i]
    )

    # Define atom ids that define vectors to move.
    moving_pairs = tuple(
        (fg.get_bromine().get_id(), fg.get_atom().get_id())
        for fg in mol.get_functional_groups()
    )
    cg = spd.Spinner(
        step_size=0.0,
        rotation_step_size=0.5,
        num_conformers=300,
        max_attempts=500,
        potential_function=FGPotential(
            moving_pairs=moving_pairs,
            target=np.array((0, 0, 1)),
        ),
    )
    conformer = cg.get_final_conformer(
        supramolecule=smolecule,
        movable_components=movable_components,
    )

    return mol.with_position_matrix(conformer.get_position_matrix())


def rotate_bbs(mol):
    smolecule = get_supramolecule(mol)
    movable_components = tuple(
        i for i, c in enumerate(smolecule.get_components())
        # Set to the known bb size.
        if c.get_num_atoms() == 20
    )

    # Define atom ids that define vectors to move.
    moving_pairs = tuple(
        tuple(i.get_id() for i in fg.get_bonders())
        for fg in mol.get_functional_groups()
    )
    cg = spd.Spinner(
        step_size=0.0,
        rotation_step_size=0.5,
        num_conformers=300,
        max_attempts=500,
        potential_function=FGPotential(
            moving_pairs=moving_pairs,
            target=np.array((0, 0, 1)),
        ),
    )
    conformer = cg.get_final_conformer(
        supramolecule=smolecule,
        movable_components=movable_components,
    )

    return mol.with_position_matrix(conformer.get_position_matrix())
