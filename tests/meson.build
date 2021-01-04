# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

d4tester = executable(
  'd4test',
  sources: files(
    'tests_dftd.f90',
    'assertion.f90',
    'class_molecule.f90',
    'geometry_reader.f90',
    'eeq_model.f90',
    'dftd4.f90',
    'pbc_tools.f90',
  ),
  dependencies: dftd4_dep_static,
)

# more specific tests are implemented by the tester binary
test('Molecule Class axis',d4tester,args: ['class_molecule','axis'])
test('Molecule Class MIC', d4tester,args: ['class_molecule','mic'])

test('Geometry Reader coord 3D',d4tester,args: ['geometry_reader','coord_3d'])
test('Geometry Reader coord 2D',d4tester,args: ['geometry_reader','coord_2d'])
test('Geometry Reader coord 1D',d4tester,args: ['geometry_reader','coord_1d'])
test('Geometry Reader coord 0D',d4tester,args: ['geometry_reader','coord_0d'])
test('Geometry Reader Xmol  0D',d4tester,args: ['geometry_reader','xmol_0d'])
test('Geometry Reader POSCAR',  d4tester,args: ['geometry_reader','poscar_3d'])

test('PBC tools convert',   d4tester,args: ['pbc_tools','convert'])
test('PBC tools cutoff',    d4tester,args: ['pbc_tools','cutoff'])

test('EEQ model water',        d4tester,args: ['eeq_model','water'])
test('EEQ model 3D Ewald',     d4tester,args: ['eeq_model','ewald'])
test('EEQ model num. gradient',d4tester,args: ['eeq_model','numgrad'])
test('EEQ model num. gradient (PBC)',d4tester,args: ['eeq_model','numgrad_pbc'])

test('Dispersion properties',    d4tester,args: ['dftd4','properties'])
test('Dispersion energies',      d4tester,args: ['dftd4','energies'])
test('Dispersion num. gradient', d4tester,args: ['dftd4','numgrad'])
test('Dispersion energies (PBC)',d4tester,args: ['dftd4','pbc_disp'])
test('Dispersion num. stress',   d4tester,args: ['dftd4','numstress'])
test('Dispersion API',           d4tester,args: ['dftd4','api'])
test('Dispersion API (PBC)',     d4tester,args: ['dftd4','pbc_api'])
test('Dispersion cell gradient', d4tester,args: ['dftd4','cell_gradient'])
test('Dispersion supercell',     d4tester,args: ['dftd4','supercell'], timeout: 120)

# currently we only test the API from the C/C++ side on Linux
if build_machine.system() == 'linux'
  if has_cpp
    test('API-C++',
      executable(
        'dftd4_cpp',
        files('cpp_api_example.cpp'),
        dependencies: dftd4_dep,
      ),
    )
  endif
  if has_cc
    test(
      'API-C',
      executable(
        'dftd4_c',
        files('c_api_example.c'),
        dependencies: dftd4_dep,
      ),
    )
  endif
endif
