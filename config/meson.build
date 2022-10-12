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
fc = meson.get_compiler('fortran')
cc = fc
if has_cc
  cc = meson.get_compiler('c')
endif
fc_id = fc.get_id()

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
endif

if get_option('openmp')
  omp_dep = dependency('openmp')
  lib_deps += omp_dep
endif

lapack_vendor = get_option('lapack')
if lapack_vendor == 'auto'
  if fc_id == 'intel'
    lapack_vendor = 'mkl'
  endif
endif

if get_option('ilp64')
  ilp64 = true
  add_project_arguments('-DIK=i8', language:'fortran')
else
  add_project_arguments('-DIK=i4', language:'fortran')
  ilp64 = false
endif

if lapack_vendor == 'mkl'
  mkl_dep = []
  if fc_id == 'intel'
    mkl_dep += cc.find_library(ilp64 ? 'mkl_intel_ilp64' : 'mkl_intel_lp64') 
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_intel_thread')
    endif
  elif fc_id == 'gcc'
    mkl_dep += cc.find_library(ilp64 ? 'mkl_gf_ilp64' : 'mkl_gf_lp64')
    if get_option('openmp')
      mkl_dep += cc.find_library('mkl_gnu_thread')
    endif
  else
    error('MKL not supported for this compiler')
  endif
  if not get_option('openmp')
    mkl_dep += cc.find_library('mkl_tbb_thread')
  endif
  mkl_dep += cc.find_library('mkl_core')
  lib_deps += mkl_dep

elif lapack_vendor == 'mkl-rt'
  mkl_dep = fc.find_library('mkl_rt')
  lib_deps += mkl_dep

elif lapack_vendor == 'openblas'
  openblas_dep = dependency( ilp64 ? 'openblas64' : 'openblas', required: false)
  if not openblas_dep.found()
    openblas_dep = fc.find_library( ilp64 ? 'openblas64' : 'openblas')
  endif
  lib_deps += openblas_dep
  if not fc.links('external dsytrs; call dsytrs(); end', dependencies: openblas_dep)
    lapack_dep = dependency(ilp64 ? 'lapack64' : 'lapack', required: false)
    if not lapack_dep.found()
      lapack_dep = fc.find_library(ilp64 ? 'lapack64' : 'lapack')
    endif
    lib_deps += lapack_dep
  endif

elif lapack_vendor == 'custom'
  custom_deps = []
  libs = get_option('custom_libraries')
  if libs[0].startswith('-L')
    foreach lib: libs
      if lib != libs[0]
        custom_deps += fc.find_library(lib, dirs: libs[0].substring(2))
      endif
    endforeach
  else
    foreach lib: libs
      custom_deps += fc.find_library(lib)
    endforeach
  endif
  if (not fc.links('external dsytrs; call dsytrs(); end', dependencies: [custom_deps,omp_dep]))
    error('Custom LAPACK libraries do not link')
  elif (not fc.links('external dsytrs; call dgemm(); end', dependencies: [custom_deps,omp_dep]))
    error('Custom BLAS libraries do not link')
  endif
  lib_deps += custom_deps

else
  lapack_dep = dependency(ilp64 ? 'lapack64' : 'lapack', required: false)
  if not lapack_dep.found()
    lapack_dep = fc.find_library(ilp64 ? 'lapack64' : 'lapack')
  endif
  lib_deps += lapack_dep
  blas_dep = fc.find_library(ilp64 ? 'blas64' : 'blas', required: false)
  if not blas_dep.found()
    blas_dep = fc.find_library(ilp64 ? 'blas64' : 'blas')
  endif
  lib_deps += blas_dep
endif

# Create the tool chain library as subproject
mctc_dep = dependency(
  'mctc-lib',
  version: '>=0.3.0',
  fallback: ['mctc-lib', 'mctc_dep'],
  default_options: ['default_library=static'],
)
lib_deps += mctc_dep

# Create the electrostatic library as subproject
multicharge_dep = dependency(
  'multicharge',
  fallback: ['multicharge', 'multicharge_dep'],
  default_options: ['default_library=static'],
)
lib_deps += multicharge_dep
