#! /usr/bin/env python3

# <<BEGIN-copyright>>
# Copyright 2022, Lawrence Livermore National Security, LLC.
# See the top-level COPYRIGHT file for details.
# 
# SPDX-License-Identifier: BSD-3-Clause
# <<END-copyright>>

description = '''
This script examines all "*.py" files in the FUDGE bin directory and prints a summary about each.
'''

import re
import pathlib
import argparse
import textwrap
import shutil

numberOfColumns = shutil.get_terminal_size().columns

parser = argparse.ArgumentParser(description=description)
args = parser.parse_args()

binDir = pathlib.Path(__file__).resolve().parent
thisScriptsName = pathlib.Path(__file__).name

scripts = sorted(binDir.glob('*.py'))

data = {}
for script in scripts:
    if script.name == thisScriptsName:
        continue
    with script.open('r') as fIn:
        lines = fIn.readlines()
        for line in lines:
            if re.match('^summaryDocStringFUDGE\s*=\s*.+$', line):
                data[script.name] = line.split('=')[-1].strip().strip("'")
                break

nameWidth = 0
for name in data:
    nameWidth = max(nameWidth, len(name))
format = '%%-%ds - %%s' % nameWidth

indent = (nameWidth + 3) * ' '

for name in data:
    line = format % (name, data[name])
    print('\n'.join(textwrap.wrap(line, width=numberOfColumns,  subsequent_indent=indent)))
