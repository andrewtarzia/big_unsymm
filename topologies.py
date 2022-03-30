#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for the m6l6 topology graph.

Author: Andrew Tarzia

"""

import math
import stk


class M6L6(stk.cage.Cage):
    """
    Represents a cage topology graph.

    Both `corner` and `linker` vertices require building blocks with
    two functional groups for this topology. This class replaces the
    `building_blocks` parameter with the `corner` and `linker`
    parameters.

    See :class:`.Cage` for more details and examples.

    """

    def __init__(
        self,
        corners,
        linkers,
        vertex_alignments=None,
        reaction_factory=stk.GenericReactionFactory(),
        num_processes=1,
        optimizer=stk.NullOptimizer(),
    ):
        """
        Initialize a :class:`.M4L4Square`.

        Parameters
        ----------
        corners : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all corner vertices on the topology
            graph.

        linkers : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all linker vertices on the topology
            graph.

        vertex_alignments : :class:`dict`, optional
            A mapping from the id of a :class:`.Vertex`
            to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is referred to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

        """

        if isinstance(corners, dict):
            building_blocks = corners
        else:
            building_blocks = {corners: (0, 1, 2, 3, 4, 5)}

        if isinstance(linkers, dict):
            linkers_dict = linkers
        else:
            linkers_dict = {linkers: (6, 7, 8, 9, 10, 11)}

        building_blocks.update(
            (building_block, vertices)
            for building_block, vertices in linkers_dict.items()
        )

        super().__init__(
            building_blocks,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
        )

    def _get_scale(self, building_block_vertices):
        return 2*super()._get_scale(building_block_vertices)

    _a = 1
    _b = math.sqrt(3)/2 * _a

    _vertex_prototypes = (
        stk.cage.LinearVertex(0, [_a, 0, 0]),
        stk.cage.LinearVertex(1, [_a/2, _b, 0]),
        stk.cage.LinearVertex(2, [-_a/2, _b, 0]),
        stk.cage.LinearVertex(3, [-_a, 0, 0]),
        stk.cage.LinearVertex(4, [-_a/2, -_b, 0]),
        stk.cage.LinearVertex(5, [_a/2, -_b, 0]),

        stk.cage.LinearVertex(6, [_b, _a/2, 0], False),
        stk.cage.LinearVertex(7, [0, _a, 0], False),
        stk.cage.LinearVertex(8, [-_b, _a/2, 0], False),
        stk.cage.LinearVertex(9, [-_b, -_a/2, 0], False),
        stk.cage.LinearVertex(10, [0, -_a, 0], False),
        stk.cage.LinearVertex(11, [_b, -_a/2, 0], False),
    )

    _edge_prototypes = (
        stk.Edge(0, _vertex_prototypes[0], _vertex_prototypes[6]),
        stk.Edge(1, _vertex_prototypes[1], _vertex_prototypes[6]),

        stk.Edge(2, _vertex_prototypes[1], _vertex_prototypes[7]),
        stk.Edge(3, _vertex_prototypes[2], _vertex_prototypes[7]),

        stk.Edge(4, _vertex_prototypes[2], _vertex_prototypes[8]),
        stk.Edge(5, _vertex_prototypes[3], _vertex_prototypes[8]),

        stk.Edge(6, _vertex_prototypes[3], _vertex_prototypes[9]),
        stk.Edge(7, _vertex_prototypes[4], _vertex_prototypes[9]),

        stk.Edge(8, _vertex_prototypes[4], _vertex_prototypes[10]),
        stk.Edge(9, _vertex_prototypes[5], _vertex_prototypes[10]),

        stk.Edge(10, _vertex_prototypes[5], _vertex_prototypes[11]),
        stk.Edge(11, _vertex_prototypes[0], _vertex_prototypes[11]),
    )

    _num_windows = 1
    _num_window_types = 1


class M12L24(stk.cage.M12L24):

    def _get_scale(self, building_block_vertices):
        return 1.5*super()._get_scale(building_block_vertices)


class OneTwoSquarePlanar(stk.metal_complex.MetalComplex):

    _metal_vertex_prototypes = (
        stk.metal_complex.MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        stk.metal_complex.MonoDentateLigandVertex(1, [2.5, 0, 0]),
        stk.metal_complex.MonoDentateLigandVertex(2, [0, 2.5, 0]),
        stk.metal_complex.BiDentateLigandVertex(3, [-2.5, -2.5, 0]),
    )

    _edge_prototypes = (
        stk.Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
            position=[1.5, 0, 0],
        ),
        stk.Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
            position=[0, 1.5, 0],
        ),
        stk.Edge(
            id=2,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=[-1.5, 0, 0],
        ),
        stk.Edge(
            id=3,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[2],
            position=[0, -1.5, 0],
        ),
    )
