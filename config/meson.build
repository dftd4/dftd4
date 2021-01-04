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

os = host_machine.system()

if os == 'windows'
  add_project_link_arguments(
    '-Wl,--allow-multiple-definition',
    language: 'fortran',
  )
endif

fc = meson.get_compiler('fortran')
fc_id = fc.get_id()
if has_cc
  has_cc = fc_id == meson.get_compiler('c').get_id()
endif
if has_cpp
  has_cpp = fc_id == meson.get_compiler('cpp').get_id()
endif

if fc_id == 'gcc'
  add_project_arguments(
    '-ffree-line-length-none',
    '-fbacktrace',
    language: 'fortran',
  )
elif fc_id == 'intel'
  add_project_arguments(
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'pgi' or fc_id == 'nvidia_hpc'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    '-traceback',
    language: 'fortran',
  )
elif fc_id == 'flang'
  add_project_arguments(
    '-Mbackslash',
    '-Mallocatable=03',
    language: 'fortran',
  )
endif


# Collect dependencies
gen_deps = []

# Collect the linear algebra backend, possible choices are MKL, OpenBLAS or netlib-BLAS
la_backend = get_option('la_backend')
if la_backend == 'mkl'
  if fc.get_id() == 'intel'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    libmkl_exe += fc.find_library('mkl_intel_thread')
  else
    libmkl_exe = [fc.find_library('mkl_gf_lp64')]
    libmkl_exe += fc.find_library('mkl_gnu_thread')
  endif
  libmkl_exe += fc.find_library('mkl_core')
  libmkl_exe += fc.find_library('iomp5')
  libmkl_sha = fc.find_library('mkl_rt')
  lib_deps += libmkl_sha
  exe_deps += libmkl_exe
elif la_backend == 'openblas'
  blas_dep = dependency('openblas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('openblas', required: true)
  endif
  gen_deps += blas_dep
  openblas_provides_lapack = fc.links(
    'external dsysv; call dsysv(); end', 
    dependencies: blas_dep,
  )
  if not openblas_provides_lapack
    lapack_dep = dependency('lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library('lapack', required: true)
    endif
    gen_deps += lapack_dep
  endif
elif la_backend == 'custom'
  foreach lib: get_option('custom_libraries')
    gen_deps += fc.find_library(lib)
  endforeach
else
  blas_dep = dependency('blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library('blas', required: true)
  endif
  gen_deps += blas_dep
  lapack_dep = dependency('lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library('lapack', required: true)
  endif
  gen_deps += lapack_dep
endif

if get_option('openmp')
  omp_dep = dependency('openmp', required: fc.get_id() != 'intel')
  if not omp_dep.found()
    omp_dep = declare_dependency(
      compile_args: '-qopenmp',
      link_args: '-qopenmp',
    )
  endif
  gen_deps += omp_dep
endif

# distribute dependencies for shared object and static executable
lib_deps += gen_deps
exe_deps += gen_deps
