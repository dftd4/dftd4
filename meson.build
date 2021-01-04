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

project(
  'dftd4',
  'fortran',
  version: '2.5.1',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')
has_cc = add_languages('c', required: false)
has_cpp = add_languages('cpp', required: false)

# General configuration information
inc_dirs = []
exe_deps = []
lib_deps = []
subdir('config')

conf = configuration_data()
conf.set('version', meson.project_version())
conf.set('commit', run_command(find_program('git'),'show','-s','--format=%h').stdout().strip())
conf.set('date', run_command(find_program('date'),'+%Y-%m-%d').stdout().strip())
conf.set('author', run_command(find_program('whoami')).stdout().strip())
conf.set('origin', run_command(find_program('hostname')).stdout().strip())
configure_file(
  input: 'source/version.f90.in',
  output: 'dftd4_version.fh',
  configuration: conf,
)

dialect = ''
bounds = ''
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  dialect = ['-std=f2008']
  bounds = '-fbounds-check'
  add_project_arguments('-fbacktrace', language: 'fortran')
endif
if fc.get_id() == 'intel'
  dialect = ['-stand','f08']
  bounds = ['-check','bounds']
  add_project_arguments('-axAVX', language: 'fortran')
  add_project_arguments('-traceback', language: 'fortran')
  add_project_link_arguments('-static', language: 'fortran')
endif

add_project_arguments(dialect, language: 'fortran')
add_project_arguments(bounds,  language: 'fortran')

# currently we only test the API from the C/C++ side on Linux
if build_machine.system() == 'linux'
  if has_cpp
    cxx = meson.get_compiler('cpp')
    if cxx.get_id() == 'gcc'
      cpp_std = '-std=c++11'
    endif
    if cxx.get_id() == 'intel'
      cpp_std = '-std=c++11'
    endif
    add_project_arguments(cpp_std, language: 'cpp')
  endif
  if has_cc
    cc = meson.get_compiler('c')
    if cc.get_id() == 'gcc'
      c_std = '-std=c11'
    endif
    if cc.get_id() == 'intel'
      c_std = '-std=c11'
    endif
    add_project_arguments(c_std, language: 'c')
  endif
endif

## ========================================== ##
## SOURCE
## ========================================== ##

fsrcs  = []
subdir('source')

## ========================================== ##
## LIBRARIES
## ========================================== ##

## ========================================== ##
## TARGETS
## ========================================== ##

# create a static library from all sources
d4lib_static = static_library(
  meson.project_name(),
  sources: fsrcs,
  include_directories: include_directories('include'),
  install: install,
  pic: true,
)

# create a shared library as well
d4lib_shared = shared_library(
  meson.project_name(),
  version: meson.project_version(),
  dependencies: lib_deps,
  link_whole: d4lib_static,
  install: install,
)

dftd4_inc = [include_directories('include'), d4lib_static.private_dir_include()]

# dependency declaration for wrap
dftd4_dep = declare_dependency(
  link_with: d4lib_shared,
  dependencies: lib_deps,
  include_directories: dftd4_inc,
)

# optional static depenendency for wrap
dftd4_dep_static = declare_dependency(
  link_with: d4lib_static,
  dependencies: exe_deps,
  include_directories: dftd4_inc,
)

# Collect main executable
subdir('app')

## ========================================== ##
## INSTALL
## ========================================== ##
if install
  install_headers('include/dftd4.h')
  
  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man1/dftd4.1.adoc'),
        output: '@BASENAME@',
      )
    )
  endif
endif

## ========================================== ##
## TESTSUITE
## ========================================== ##
subdir('tests')
subdir('python')
