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

dftd4_exe = executable(
  meson.project_name(),
  sources: files('main.f90'),
  dependencies: dftd4_dep,
  install: install,
)

test('app-version', dftd4_exe, args: '--version')
test('app-help', dftd4_exe, args: '--help')
test('app-license', dftd4_exe, args: '--license')
test('app-citation', dftd4_exe, args: '--citation')
test('app-noargs', dftd4_exe, should_fail: true)

test(
  'app-energy',
  files('tester.py'),
  args: [
    dftd4_exe,
    files('energy.json'),
    files('ammonia.tmol'),
    '--noedisp',
    '--func', 'tpssh'
  ],
)
test(
  'app-gradient',
  files('tester.py'),
  args: [
    dftd4_exe,
    files('gradient.json'),
    '--noedisp',
    '--func', 'scan',
    '--grad',
    '--',
    files('nitralin.mol'),
  ],
)
test(
  'app-pair-analysis',
  files('tester.py'),
  args: [
    dftd4_exe,
    files('pair-analysis.json'),
    files('caffeine.xyz'),
    '--noedisp',
    '--pair-resolved',
    '--func', 'b3lyp',
  ],
)
