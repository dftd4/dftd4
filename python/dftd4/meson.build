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

cc = meson.get_compiler('c')

pymod = import('python')
python = pymod.find_installation(
  'python@0@'.format(get_option('python_version')),
  modules: [
    'cffi',
  ],
)
python_dep = python.dependency(required: true)

# Python's CFFI is horrible in working with preprocessor statements,
# therefore, we have to preprocess the header before passing it to the ffibuilder
dftd4_pp = configure_file(
  command: [cc, '-E', '@INPUT@'],
  input: dftd4_header,
  output: '_libdftd4.h',
  capture: true,
)

# This is the actual out-of-line API processing of the ffibuilder
dftd4_cffi_srcs = configure_file(
  command: [python, files('..'/'ffibuilder.py'), '@INPUT@', '@BASENAME@'],
  input: dftd4_pp,
  output: '@BASENAME@.c',
)

# Actual generation of the Python extension, since the shared_module does not work
# well with dependency objects, we will trick it by linking a whole static lib
dftd4_pyext = python.extension_module(
  '_libdftd4',
  link_whole: static_library(
    '_libdftd4',
    dftd4_cffi_srcs,
    dependencies: [dftd4_dep, python_dep],
  ),
  dependencies: [dftd4_dep, python_dep],
  install: install,
  subdir: 'dftd4',
)

if install
  python.install_sources(
    files(
      '__init__.py',
      'ase.py',
      'interface.py',
      'libdftd4.py',
      'qcschema.py',
      'test_ase.py',
      'test_interface.py',
      'test_libdftd4.py',
      'test_qcschema.py',
    ),
    subdir: 'dftd4',
  )
endif

