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

mypy = find_program('mypy', required: false)
pytest = find_program('pytest', required: false)
pylint = find_program('pylint', required: false)

py_dftd4 = meson.current_source_dir() / 'dftd4'
pypath = environment()
pypath.prepend('LD_LIBRARY_PATH', meson.current_build_dir() / '..')
pypath.prepend('PYTHONPATH', meson.current_source_dir())

if pytest.found()
  test('pytest: dftd4.py', pytest, args: ['--pyargs', 'dftd4'], env: pypath)
endif

if mypy.found()
  test('mypy: dftd4.py', mypy, args: [py_dftd4, '--ignore-missing-imports'],
       env: pypath)
endif

if pylint.found()
  test('pylint: dftd4.py', pylint, args: [py_dftd4, '--disable=invalid-names'],
       env: pypath)
endif
