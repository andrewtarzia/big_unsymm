#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for the cis-protected square planar topology graph.

Author: Andrew Tarzia

"""

import stk


class CisProtectedMonodentate(stk.metal_complex.MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    """

    _metal_vertex_prototypes = (
        stk.metal_complex.MetalVertex(0, [0, 0, 0]),
    )
    _ligand_vertex_prototypes = (
        stk.metal_complex.MonoDentateLigandVertex(1, [2.5, 0, 0]),
        stk.metal_complex.MonoDentateLigandVertex(2, [0, 2.5, 0]),
    )

    _edge_prototypes = (
        stk.Edge(
            id=0,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[0],
        ),
        stk.Edge(
            id=1,
            vertex1=_metal_vertex_prototypes[0],
            vertex2=_ligand_vertex_prototypes[1],
        ),
    )
