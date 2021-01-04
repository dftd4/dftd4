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
  dependencies: dftd4_dep_static,
  install: not meson.is_subproject(),
)

# some very basic checks to see if the executable reacts at all
test('Print version', dftd4_exe, args: '--version')
test('Print help', dftd4_exe, args: '--help')
test('Print license', dftd4_exe, args: '--license')
test('No arguments', dftd4_exe, should_fail: true)
