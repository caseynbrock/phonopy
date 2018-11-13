# Copyright (C) 2015 Atsushi Togo
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import sys
import numpy as np

from phonopy.file_IO import collect_forces
from phonopy.interface.vasp import (get_scaled_positions_lines,
                                    sort_positions_by_symbols,
                                    check_forces,
                                    get_drift_forces)
from phonopy.units import Bohr
from phonopy.structure.atoms import PhonopyAtoms as Atoms
from phonopy.structure.atoms import symbol_map


def parse_set_of_forces(num_atoms, forces_filenames, verbose=True):
    hook = 'Atomic forces:'
    is_parsed = True
    force_sets = []

    for i, filename in enumerate(forces_filenames):
        if verbose:
            sys.stdout.write("%d. " % (i + 1))
        with open(filename, 'r') as fin:
            f = fin.readlines()
        socorro_forces = collect_forces(f,
                                    num_atoms,
                                    hook,
                                    [1,2,3],
                                    skiplines=2)
        if check_forces(socorro_forces, num_atoms, filename, verbose=verbose):
            drift_force = get_drift_forces(socorro_forces,
                                           filename=filename,
                                           verbose=verbose)
            force_sets.append(np.array(socorro_forces) - drift_force)
        else:
            is_parsed = False

    if is_parsed:
        return force_sets
    else:
        return []

def read_socorro(filename):
    """
    Read file with crystal structure in socorro format. 
    Return PhonopyAtoms instance.
    """
    socorro_in = SocorroIn(open(filename).readlines())
    tags = socorro_in.get_variables()
    avec = [tags['scale'][i] * np.array(tags['avec'][i]) for i in range(3)]

    symbols = tags['atoms']['spfnames']
    numbers = []
    for s in symbols:
        numbers.append(symbol_map[s])
            
    return Atoms(numbers=numbers,
                 cell=avec,
                 scaled_positions=tags['atoms']['positions'])

def write_socorro(filename, cell):
    with open(filename, 'w') as f:
        f.write(get_socorro_structure(cell))

def write_supercells_with_displacements(supercell,
                                        cells_with_displacements):
    write_socorro("supercell.in", supercell)
    for i, cell in enumerate(cells_with_displacements):
        write_socorro("supercell-%03d.in" % (i + 1), cell)

def get_socorro_structure(cell):
    lattice = cell.get_cell()
    positions = cell.get_scaled_positions()
    chemical_symbols = cell.get_chemical_symbols()

    lines = "title\n"
    lines += "1.0\n"
    lines += ((" %21.16f" * 3 + "\n") * 3) % tuple(lattice.ravel())
    lines += "lattice\n"
    lines += "%d\n" % len(chemical_symbols)
    for pos, symbol in zip(positions,chemical_symbols):
        lines += ("%s " + "%21.16lf"*3 +"\n") % tuple([symbol] + pos.tolist())
        
    return lines

class SocorroIn(object):
    """
    Inputs:
        lines: socorro formatted crystal file containing structural parameters
    
    Public methods:
        get_variables(): returns a library of structure parameters read in
    """
    def __init__(self, lines):
        self._set_methods = {'scale':  self._set_scale,
                             'avec':   self._set_avec,
                             'atoms':  self._set_atoms}
        self._tags = {'atoms': None,
                      'avec':  None,
                      'scale': [1.0, 1.0, 1.0]}
        self._lines = iter(lines)
        self._collect()

    def get_variables(self):
        return self._tags

    def _collect(self):
        """
        socorro has a rigid crystal file structure that does not require advanced parsing
        """
        self._lines.next()
        self._set_scale()
        self._set_avec()
        self._lines.next()
        self._set_atoms()

    def _set_atoms(self):
        natoms = int(self._lines.next().split()[0])
        spfnames = []
        positions = []
        for j in range(natoms):            
            line = self._lines.next()
            spfnames.append(line.split()[0])
            positions.append(
                [float(x) for x in line.split()[1:4]])
            
        self._tags['atoms'] = {'spfnames':  spfnames,
                               'positions': positions}
            
    def _set_avec(self):
        avec = []
        for i in range(3):
            avec.append([float(x) for x in self._lines.next().split()[:3]])
        self._tags['avec'] = avec
            
        self._tags['atoms'] = {'spfnames':  spfnames,
                               'positions': positions}
            
    def _set_avec(self):
        avec = []
        for i in range(3):
            avec.append([float(x) for x in self._lines.next().split()[:3]])
        self._tags['avec'] = avec

    def _set_scale(self):
        scale = float(self._lines.next().split()[0])
        for i in range(3):
            self._tags['scale'][i] = scale

if __name__ == '__main__':
    import sys
    from phonopy.structure.symmetry import Symmetry
    cell, sp_filenames = read_socorro(sys.argv[1])
    symmetry = Symmetry(cell)
    print("# %s" % symmetry.get_international_table())
    print(get_socorro_structure(cell, sp_filenames=sp_filenames))
