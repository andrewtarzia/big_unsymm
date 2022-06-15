#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to build all cages in this project.

Author: Andrew Tarzia

"""

import logging
import sys
import os
import glob
import shutil
import subprocess as sp

from env_set import pd24_path, calc_path, xtb_path


def pdb_to_xyz(pdb, xyz):

    xyz_string = [None, 'been cleaning\n']

    with open(pdb, 'r') as f:
        lines = f.readlines()

    atom_count = 0
    for line in lines:
        if 'HETATM' in line:
            l = line.rstrip().split()
            atype = l[2]
            x = float(l[5])
            y = float(l[6])
            z = float(l[7])
            xyz_line = f'{atype} {x} {y} {z}\n'
            atom_count += 1
            xyz_string.append(xyz_line)

    xyz_string[0] = f'{atom_count}\n'
    with open(xyz, 'w') as f:
        f.write(''.join(xyz_string))

def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 0 arguments:'
        )
        sys.exit()
    else:
        pass

    _wd = pd24_path()
    _cd = calc_path()
    init_dir = os.getcwd()

    if not os.path.exists(_wd):
        os.mkdir(_wd)

    if not os.path.exists(_cd):
        os.mkdir(_cd)

    unopt_files = glob.glob(str(_wd / '*.pdb'))
    for unopt_file in unopt_files:
        unopt_pdb = unopt_file.split('/')[-1]
        name = unopt_pdb.replace('.pdb', '')
        cleaned_file = _wd / f'{name}_cleaned.xyz'
        opt_file = _wd / f'{name}_opt.xyz'
        charge = 48
        pdb_to_xyz(unopt_file, cleaned_file)

        output_dir = _wd / f'{name}_xtbopt'
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)
        os.mkdir(output_dir)
        os.chdir(output_dir)
        if not os.path.exists(opt_file):
            try:
                out_file = 'optimisation.out'
                logging.info(f'optimising {name}')
                # Modify the memory limit.
                memory = 'ulimit -s unlimited ;'
                # Set optimization level and type.
                optimization = f'--opt normal'

                cmd = (
                    f'{memory} {xtb_path()} {cleaned_file} '
                    f'--gfn 2 '
                    f'{optimization} --parallel 4 '
                    f'--chrg {charge} '
                )
                logging.info(f'running {cmd}')

                with open(out_file, 'w') as f:
                    # Note that sp.call will hold the program until completion
                    # of the calculation.
                    sp.call(
                        cmd,
                        stdin=sp.PIPE,
                        stdout=f,
                        stderr=sp.PIPE,
                        # Shell is required to run complex arguments.
                        shell=True
                    )

                os.system(f'cp xtbopt.xyz {opt_file}')
            finally:
                os.chdir(init_dir)

        with open(opt_file, 'r') as f:
            line = f.readlines()[1]
        print(line)

        energy = 0
        raise SystemExit('get ey')
        logging.info(f'{name} optimised energy: {energy}')


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.DEBUG,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
