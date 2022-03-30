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
        if 'tri' in struct:
            continue
        if 'sqr' in struct:
            continue
        if 'hex' in struct:
            continue
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
    ax.set_xlabel('ligand', fontsize=16)
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
    ax.set_xlabel('ligand', fontsize=16)
    if per_ligand:
        ax.set_ylabel(r'rel. energy [kJ mol$^{-1}$]', fontsize=16)
    else:
        ax.set_ylabel(
            r'rel. energy per L [kJ mol$^{-1}$]', fontsize=16
        )
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
