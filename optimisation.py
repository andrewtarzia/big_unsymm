#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for optimisation functions.

Author: Andrew Tarzia

"""

import logging
import os
import stko

import env_set


def MOC_MD_opt(
    cage,
    cage_name,
    integrator,
    temperature,
    N,
    timestep,
    equib,
    production,
    opt_conf,
    metal_FFs,
    metal_ligand_bond_order='half',
    save_conf=False,
    gulp_exec=None,
):
    """
    Perform UFF4MOF molecular dynamics of MOC.
    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.
    cage_name : :class:`str`
        Name of cage.
    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.
    """
    raise NotImplementedError('no print')

    if gulp_exec is None:
        gulp_exec = '/home/atarzia/software/gulp-5.1/Src/gulp/gulp'

    print(f'..........doing UFF4MOF MD of {cage_name}')
    gulp_MD = stko.GulpUFFMDOptimizer(
        gulp_path=gulp_exec,
        metal_FF=metal_FFs,
        metal_ligand_bond_order=metal_ligand_bond_order,
        output_dir=f'cage_opt_{cage_name}_MD',
        integrator=integrator,
        ensemble='nvt',
        temperature=temperature,
        equilbration=equib,
        production=production,
        timestep=timestep,
        N_conformers=N,
        opt_conformers=opt_conf,
        save_conformers=save_conf
    )
    gulp_MD.assign_FF(cage)
    cage = gulp_MD.optimize(cage)

    return cage


def MOC_xtb_conformers(
    cage,
    cage_name,
    etemp,
    output_dir,
    conformer_dir,
    nc,
    free_e,
    charge,
    gfn_exec=None,
    opt=False,
    opt_level=None,
    solvent=None,
    handle_failure=False
):
    """
    Perform GFN2-xTB conformer scan of MOC.
    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.
    cage_name : :class:`str`
        Name of cage.
    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.
    """
    raise NotImplementedError('no print')

    if gfn_exec is None:
        gfn_exec = '/home/atarzia/software/xtb-190806/bin/xtb'

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print(
        f'..........doing XTB conformer sorting by '
        f'energy of {cage_name}'
    )
    conformers = glob.glob(f'{conformer_dir}/conf_*.xyz')
    ids = []
    energies = []
    min_energy = 10E20
    for file in sorted(conformers):
        id = file.replace('.xyz', '').split('_')[-1]
        cage = cage.with_structure_from_file(file)
        opt_failed = False
        if opt:
            print(f'optimising conformer {id}')
            xtb_opt = stko.XTB(
                xtb_path=gfn_exec,
                output_dir=f'opt_{cage_name}_{id}',
                gfn_version=2,
                num_cores=nc,
                opt_level=opt_level,
                charge=charge,
                num_unpaired_electrons=free_e,
                max_runs=1,
                electronic_temperature=etemp,
                calculate_hessian=False,
                unlimited_memory=True,
                solvent=solvent_str,
                solvent_grid=solvent_grid
            )
            try:
                cage = xtb_opt.optimize(mol=cage)
                cage.write(os.path.join(
                    f'{output_dir}',
                    f'conf_{id}_opt.xyz',
                ))
            except stko.XTBConvergenceError:
                if handle_failure:
                    opt_failed = True
                    print(f'optimising conformer {id}: FAILED')
                else:
                    raise stko.XTBConvergenceError()

        print(f'..........calculating energy of {id} of {cage_name}')
        # Extract energy.
        xtb_energy = stko.XTBEnergy(
            xtb_path=gfn_exec,
            output_dir=f'ey_{cage_name}_{id}',
            num_cores=nc,
            charge=charge,
            num_unpaired_electrons=free_e,
            electronic_temperature=etemp,
            unlimited_memory=True,
            solvent=solvent_str,
            solvent_grid=solvent_grid
        )
        if handle_failure and opt_failed:
            energy = 10E24
        else:
            energy = xtb_energy.get_energy(cage)
        if energy < min_energy:
            min_energy_conformer = file
            min_energy = energy
        ids.append(id)
        energies.append(energy)

    print('done', min_energy, min_energy_conformer)
    cage = cage.with_structure_from_file(min_energy_conformer)

    energies = [(i-min(energies))*2625.5 for i in energies]
    fig, ax = scatter_plot(
        X=ids, Y=energies,
        xtitle='conformer id',
        ytitle='rel. energy [kJmol$^{-1}$]',
        xlim=(0, 201),
        ylim=(-5, 1000)
    )

    fig.tight_layout()
    fig.savefig(
        os.path.join(output_dir, f'{cage_name}_conf_energies.pdf'),
        dpi=720,
        bbox_inches='tight'
    )
    plt.close()

    return cage


def MOC_xtb_opt(
    cage,
    cage_name,
    nc,
    opt_level,
    etemp,
    charge,
    free_e,
    gfn_exec=None,
    solvent=None
):
    """
    Perform GFN2-xTB optimisation of MOC.
    Parameters
    ----------
    cage : :class:`stk.ConstructedMolecule`
        Cage to be optimised.
    cage_name : :class:`str`
        Name of cage.
    Returns
    -------
    cage : :class:`stk.ConstructedMolecule`
        Optimised cage.
    """
    raise NotImplementedError('no print')

    if gfn_exec is None:
        gfn_exec = '/home/atarzia/software/xtb-190806/bin/xtb'

    if solvent is None:
        solvent_str = None
        solvent_grid = 'normal'
    else:
        solvent_str, solvent_grid = solvent

    print(f'..........doing XTB optimisation of {cage_name}')
    xtb_opt = stko.XTB(
        xtb_path=gfn_exec,
        output_dir=f'cage_opt_{cage_name}_xtb',
        gfn_version=2,
        num_cores=nc,
        opt_level=opt_level,
        charge=charge,
        num_unpaired_electrons=free_e,
        max_runs=1,
        electronic_temperature=etemp,
        calculate_hessian=False,
        unlimited_memory=True,
        solvent=solvent_str,
        solvent_grid=solvent_grid
    )
    cage = xtb_opt.optimize(mol=cage)

    return cage


def optimisation_sequence(mol, name, charge, calc_dir):
    gulp1_output = os.path.join(calc_dir, f'{name}_gulp1.mol')
    gulp2_output = os.path.join(calc_dir, f'{name}_gulp2.mol')

    if not os.path.exists(gulp1_output):
        output_dir = os.path.join(calc_dir, f'{name}_gulp1')
        CG = True
        logging.info(f'UFF4MOF optimisation 1 of {name} CG: {CG}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=env_set.gulp_path(),
            maxcyc=1000,
            metal_FF={46: 'Pd4+2'},
            metal_ligand_bond_order='',
            output_dir=output_dir,
            conjugate_gradient=CG,
        )
        gulp_opt.assign_FF(mol)
        gulp1_mol = gulp_opt.optimize(mol=mol)
        gulp1_mol.write(gulp1_output)
    else:
        gulp1_mol = mol.with_structure_from_file(gulp1_output)

    if not os.path.exists(gulp2_output):
        output_dir = os.path.join(calc_dir, f'{name}_gulp2')
        CG = False
        logging.info(f'UFF4MOF optimisation 2 of {name} CG: {CG}')
        gulp_opt = stko.GulpUFFOptimizer(
            gulp_path=env_set.gulp_path(),
            maxcyc=1000,
            metal_FF={46: 'Pd4+2'},
            metal_ligand_bond_order='',
            output_dir=output_dir,
            conjugate_gradient=CG,
        )
        gulp_opt.assign_FF(gulp1_mol)
        gulp2_mol = gulp_opt.optimize(mol=gulp1_mol)
        gulp2_mol.write(gulp2_output)
    else:
        gulp2_mol = mol.with_structure_from_file(gulp2_output)

    logging.warning('update this:')
    final_mol = mol.with_structure_from_file(gulp2_output)
    return final_mol