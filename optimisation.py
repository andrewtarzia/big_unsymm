#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Module for optimisation functions.

Author: Andrew Tarzia

"""

import logging
import os
import stk
import stko
import glob

import env_set


def xtb_conformer_opt(
    mol,
    output_dir,
    conformer_dir,
    charge,
):
    """
    Perform GFN2-xTB conformer scan of MOC.

    """

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    conformers = glob.glob(f'{conformer_dir}/conf_*.xyz')

    min_energy = 10E20
    for c_file in sorted(conformers):
        id_ = c_file.replace('.xyz', '').split('_')[-1]
        conf = mol.with_structure_from_file(c_file)

        opt_conf_file = os.path.join(output_dir, f'conf_{id_}_opt.xyz')
        if not os.path.exists(opt_conf_file):
            logging.info(f'optimising conformer ID {id_} from {c_file}')
            xtb_opt = stko.XTB(
                xtb_path=env_set.xtb_path(),
                output_dir=os.path.join(output_dir, f'conf_{id_}_opt'),
                gfn_version=2,
                num_cores=4,
                opt_level='crude',
                charge=charge,
                max_runs=1,
                calculate_hessian=False,
                unlimited_memory=True,
            )
            opt_conf = xtb_opt.optimize(mol=conf)
            opt_conf.write(opt_conf_file)
        else:
            opt_conf = conf.with_structure_from_file(opt_conf_file)

        logging.info(f'calculating energy of {id_}')
        # Extract energy.
        xtb_energy = stko.XTBEnergy(
            xtb_path=env_set.xtb_path(),
            output_dir=os.path.join(output_dir, f'conf_{id_}_ey'),
            num_cores=4,
            charge=charge,
            unlimited_memory=True,
        )
        energy = xtb_energy.get_energy(opt_conf)
        if energy < min_energy:
            min_energy_conformer = (
                stk.BuildingBlock.init_from_molecule(opt_conf)
            )
            min_energy = energy

    logging.info(
        f'lowest energy conformer ID {id_} with ey: {min_energy} a.u.'
    )

    return min_energy_conformer


def optimisation_sequence(mol, name, charge, calc_dir):
    gulp1_output = os.path.join(calc_dir, f'{name}_gulp1.mol')
    gulp2_output = os.path.join(calc_dir, f'{name}_gulp2.mol')
    gulpmd_output = os.path.join(calc_dir, f'{name}_gulpmd.mol')
    xtbconfs_output = os.path.join(calc_dir, f'{name}_xtbconf.mol')
    xtbopt_output = os.path.join(calc_dir, f'{name}_xtb.mol')

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

    if not os.path.exists(gulpmd_output):
        output_dir = os.path.join(calc_dir, f'{name}_gulpmd')

        logging.info(f'UFF4MOF equilib MD of {name}')
        gulp_MD = stko.GulpUFFMDOptimizer(
            gulp_path=env_set.gulp_path(),
            metal_FF={46: 'Pd4+2'},
            metal_ligand_bond_order='',
            output_dir = os.path.join(calc_dir, f'{name}_gulpmd'),
            integrator='leapfrog verlet',
            ensemble='nvt',
            temperature=500,
            timestep=0.25,
            equilbration=0.5,
            production=0.5,
            N_conformers=2,
            opt_conformers=False,
            save_conformers=False,
        )
        gulp_MD.assign_FF(gulp2_mol)
        gulpmd_mol = gulp_MD.optimize(mol=gulp2_mol)

        logging.info(f'UFF4MOF production MD of {name}')
        gulp_MD = stko.GulpUFFMDOptimizer(
            gulp_path=env_set.gulp_path(),
            metal_FF={46: 'Pd4+2'},
            metal_ligand_bond_order='',
            output_dir = os.path.join(calc_dir, f'{name}_gulpmd'),
            integrator='leapfrog verlet',
            ensemble='nvt',
            temperature=500,
            timestep=0.75,
            equilbration=0.5,
            production=30.0,
            N_conformers=30,
            opt_conformers=True,
            save_conformers=False,
        )
        gulp_MD.assign_FF(gulpmd_mol)
        gulpmd_mol = gulp_MD.optimize(mol=gulpmd_mol)
        gulpmd_mol.write(gulpmd_output)
    else:
        gulpmd_mol = mol.with_structure_from_file(gulpmd_output)

    # if not os.path.exists(xtbconfs_output):
    #     output_dir = os.path.join(calc_dir, f'{name}_xtbconfs')
    #     conformer_dir = os.path.join(calc_dir, f'{name}_gulpmd')
    #     logging.info(f'xtb conformer ranking of {name}')
    #     xtb_conf_mol = xtb_conformer_opt(
    #         mol=gulpmd_mol,
    #         output_dir=output_dir,
    #         conformer_dir=conformer_dir,
    #         charge=charge,
    #     )
    #     xtb_conf_mol.write(xtbconfs_output)
    # else:
    #     xtb_conf_mol = mol.with_structure_from_file(xtbconfs_output)

    if not os.path.exists(xtbopt_output):
        output_dir = os.path.join(calc_dir, f'{name}_xtbopt')
        logging.info(f'final xtb optimisation of {name}')
        xtb_opt = stko.XTB(
            xtb_path=env_set.xtb_path(),
            output_dir=output_dir,
            gfn_version=2,
            num_cores=6,
            charge=charge,
            opt_level='extreme',
            num_unpaired_electrons=0,
            max_runs=1,
            calculate_hessian=False,
            unlimited_memory=True,
            solvent=None,
        )
        xtbopt_mol = xtb_opt.optimize(mol=gulpmd_mol)
        xtbopt_mol.write(xtbopt_output)
    else:
        xtbopt_mol = mol.with_structure_from_file(xtbopt_output)

    final_mol = mol.with_structure_from_file(xtbopt_output)
    return final_mol