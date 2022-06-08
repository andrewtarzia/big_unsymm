#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for plotting functions.

Author: Andrew Tarzia

"""

import matplotlib.pyplot as plt
import os

from env_set import figu_path


def plot_strain_energies(results_dict, outname):
    fig, ax = plt.subplots(figsize=(8, 5))

    x_position = 0
    _x_names = []
    all_values = []
    for struct in results_dict:
        s_values = results_dict[struct]
        x_position += 1
        all_values.append((x_position, s_values['sum_strain_energy']))
        _x_names.append((x_position, struct))

    ax.scatter(
        x=[i[0] for i in all_values],
        y=[i[1] for i in all_values],
        c='gold',
        edgecolors='k',
        s=180,
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('structure', fontsize=16)
    ax.set_ylabel(r'sum strain energy [kJ mol$^{-1}$]', fontsize=16)
    # ax.set_xlim((0, 1))
    ax.set_ylim(0, None)
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])

    fig.tight_layout()
    fig.savefig(
        os.path.join(figu_path(), f'{outname}.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_ops(results_dict, outname):

    fig, ax = plt.subplots(figsize=(8, 5))

    x_position = 0
    _x_names = []
    all_values = []
    for struct in results_dict:
        s_values = results_dict[struct]
        x_position += 1
        all_values.append((x_position, s_values['min_order_param']))
        _x_names.append((x_position, struct))

    ax.scatter(
        x=[i[0] for i in all_values],
        y=[i[1] for i in all_values],
        c='gold',
        edgecolors='k',
        s=180,
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('ligand', fontsize=16)
    ax.set_ylabel(r'min op', fontsize=16)
    # ax.set_xlim((0, 1))
    ax.set_ylim(0, 1)
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])

    fig.tight_layout()
    fig.savefig(
        os.path.join(figu_path(), f'{outname}.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_energies(results_dict, outname, per_ligand=False):

    fig, ax = plt.subplots(figsize=(8, 5))

    x_position = 0
    _x_names = []
    all_energies = []
    for struct in results_dict:
        s_values = results_dict[struct]
        if 'tri' in struct:
            y = s_values['xtb_energy']/3
        elif 'sqr' in struct:
            y = s_values['xtb_energy']/4
        elif 'hex' in struct:
            y = s_values['xtb_energy']/6
        else:
            y = s_values['xtb_energy']

        x_position += 1
        all_energies.append((x_position, y))
        _x_names.append((x_position, struct))

    min_energy = min([i[1] for i in all_energies])
    print(all_energies, min_energy)
    ax.scatter(
        x=[i[0] for i in all_energies],
        y=[(i[1]-min_energy)*2625.5 for i in all_energies],
        c='gold',
        edgecolors='k',
        s=180,
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('structure', fontsize=16)
    if per_ligand:
        ax.set_ylabel(
            r'rel. energy per L [kJ mol$^{-1}$]', fontsize=16
        )
    else:
        ax.set_ylabel(r'rel. energy [kJ mol$^{-1}$]', fontsize=16)

    # ax.set_xlim((0, 1))
    ax.set_ylim(-0.1, None)
    ax.set_xticks([i[0] for i in _x_names])
    ax.set_xticklabels([i[1] for i in _x_names])

    fig.tight_layout()
    fig.savefig(
        os.path.join(figu_path(), f'{outname}.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_subs_energies(results_dict, outname):

    fig, ax = plt.subplots(figsize=(8, 5))

    tri_energies = []
    sqr_energies = []
    hex_energies = []
    for struct in results_dict:
        s_values = results_dict[struct]
        y = s_values['xtb_energy']
        x_position = int(struct[-1])
        if 'tri' in struct:
            tri_energies.append((x_position, y))
        elif 'sqr' in struct:
            sqr_energies.append((x_position, y))
        elif 'hex' in struct:
            hex_energies.append((x_position, y))

    min_tri_energy = min([i[1] for i in tri_energies])
    min_sqr_energy = min([i[1] for i in sqr_energies])
    min_hex_energy = min([i[1] for i in hex_energies])
    print(tri_energies, min_tri_energy)
    print(sqr_energies, min_sqr_energy)
    print(hex_energies, min_hex_energy)
    ax.scatter(
        x=[i[0] for i in tri_energies],
        y=[(i[1]-min_tri_energy)*2625.5 for i in tri_energies],
        c='gold',
        edgecolors='k',
        s=180,
        label='tri'
    )
    ax.scatter(
        x=[i[0] for i in sqr_energies],
        y=[(i[1]-min_sqr_energy)*2625.5 for i in sqr_energies],
        c='skyblue',
        edgecolors='k',
        s=180,
        label='sqr'
    )
    ax.scatter(
        x=[i[0] for i in hex_energies],
        y=[(i[1]-min_hex_energy)*2625.5 for i in hex_energies],
        c='forestgreen',
        edgecolors='k',
        s=180,
        label='hex'
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('structure id', fontsize=16)
    ax.set_ylabel(r'rel. energy [kJ mol$^{-1}$]', fontsize=16)

    # ax.set_xlim((0, 1))
    ax.set_ylim(-0.1, None)
    # ax.set_xticks([i[0] for i in _x_names])
    # ax.set_xticklabels([i[1] for i in _x_names])
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(figu_path(), f'{outname}.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()


def plot_man_subs_energies(results_dict, outname):

    fig, ax = plt.subplots(figsize=(8, 5))

    A1_energies = []
    A2_energies = []
    A3_energies = []
    for struct in results_dict:
        y = results_dict[struct]
        x_position = int(struct[-1])
        if 'A1' in struct:
            A1_energies.append((x_position, y))
        elif 'A2' in struct:
            A2_energies.append((x_position, y))
        elif 'A3' in struct:
            A3_energies.append((x_position, y))

    min_energy = min([
        min([i[1] for i in A1_energies]),
        min([i[1] for i in A2_energies]),
        min([i[1] for i in A3_energies]),
    ])
    print(A1_energies, min_energy)
    print(A2_energies)
    print(A3_energies)
    ax.scatter(
        x=[i[0] for i in A1_energies],
        y=[(i[1]-min_energy)*2625.5 for i in A1_energies],
        c='gold',
        edgecolors='k',
        s=180,
        label='A1'
    )
    ax.scatter(
        x=[i[0] for i in A2_energies],
        y=[(i[1]-min_energy)*2625.5 for i in A2_energies],
        c='skyblue',
        edgecolors='k',
        s=180,
        label='A2'
    )
    ax.scatter(
        x=[i[0] for i in A3_energies],
        y=[(i[1]-min_energy)*2625.5 for i in A3_energies],
        c='forestgreen',
        edgecolors='k',
        s=180,
        label='A3'
    )

    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.set_xlabel('structure id', fontsize=16)
    ax.set_ylabel(r'rel. energy [kJ mol$^{-1}$]', fontsize=16)

    # ax.set_xlim((0, 1))
    ax.set_ylim(-0.1, None)
    # ax.set_xticks([i[0] for i in _x_names])
    # ax.set_xticklabels([i[1] for i in _x_names])
    ax.legend(fontsize=16)

    fig.tight_layout()
    fig.savefig(
        os.path.join(figu_path(), f'{outname}.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()
