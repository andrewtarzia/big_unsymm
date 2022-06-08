#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Distributed under the terms of the MIT License.

"""
Script to analyse all subsystems extracted.

Author: Andrew Tarzia

"""

import logging
import glob
import sys
import os
import shutil
import re
import subprocess as sp

from env_set import manu_subs_path, xtb_path
from plotting import plot_man_subs_energies


def main():
    if (not len(sys.argv) == 1):
        logging.info(
            f'Usage: {__file__}\n'
            '   Expected 1 arguments:'
        )
        sys.exit()
    else:
        pass

    _wd = manu_subs_path()

    squares = sorted(glob.glob(os.path.join(_wd, '*_sqr_*.xyz')))
    triangles = sorted(glob.glob(os.path.join(_wd, '*_tri_*.xyz')))

    logging.info(
        f'there are {len(squares) + len(triangles)} structures.'
    )

    square_results = {}
    triangle_results = {}
    for s_file in squares:
        # Make a dir.
        name = s_file.split('/')[-1].replace('.xyz', '')
        xyz = s_file.split('/')[-1]
        dire = _wd / f'{name}_xtbey'
        xtb_output_file = dire / 'energy.output'

        if not os.path.exists(xtb_output_file):
            if os.path.exists(dire):
                shutil.rmtree(dire)
            os.mkdir(dire)
            os.chdir(dire)
            os.system(f'cp {_wd / xyz} .')
            try:
                # Run xtb command.
                charge = 8
                cmd = (
                    f'ulimit -s unlimited ; {xtb_path()} {xyz} '
                    f'--parallel 4 --chrg {charge}'
                )
                logging.info(f'running {cmd}')
                with open('energy.output', 'w') as f:
                    # Note that sp.call will hold the program until
                    # completion of the calculation.
                    sp.call(
                        cmd,
                        stdin=sp.PIPE,
                        stdout=f,
                        stderr=sp.PIPE,
                        # Shell is required to run complex arguments.
                        shell=True
                    )
            finally:
                os.chdir(_wd)

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")

        with open(xtb_output_file, 'r') as f:
            for line in f.readlines():
                if '          | TOTAL ENERGY  ' in line:
                    string = nums.search(line.rstrip()).group(0)
                    energy = float(string)
                    break

        square_results[name] = energy

    for s_file in triangles:
        # Make a dir.
        name = s_file.split('/')[-1].replace('.xyz', '')
        xyz = s_file.split('/')[-1]
        dire = _wd / f'{name}_xtbey'
        xtb_output_file = dire / 'energy.output'

        if not os.path.exists(xtb_output_file):
            if os.path.exists(dire):
                shutil.rmtree(dire)
            os.mkdir(dire)
            os.chdir(dire)
            os.system(f'cp {_wd / xyz} .')
            try:
                # Run xtb command.
                charge = 6
                cmd = (
                    f'ulimit -s unlimited ; {xtb_path()} {xyz} '
                    f'--parallel 4 --chrg {charge}'
                )
                logging.info(f'running {cmd}')
                with open('energy.output', 'w') as f:
                    # Note that sp.call will hold the program until
                    # completion of the calculation.
                    sp.call(
                        cmd,
                        stdin=sp.PIPE,
                        stdout=f,
                        stderr=sp.PIPE,
                        # Shell is required to run complex arguments.
                        shell=True
                    )
            finally:
                os.chdir(_wd)

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")

        with open(xtb_output_file, 'r') as f:
            for line in f.readlines():
                if '          | TOTAL ENERGY  ' in line:
                    string = nums.search(line.rstrip()).group(0)
                    energy = float(string)
                    break

        triangle_results[name] = energy

    plot_man_subs_energies(
        results_dict=square_results,
        outname='man_sqr_subs_energies',
    )
    plot_man_subs_energies(
        results_dict=triangle_results,
        outname='man_tri_subs_energies',
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s | %(levelname)s | %(message)s',
    )
    main()
