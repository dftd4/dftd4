# 
project('dftd4', 'fortran', 'cpp',
        version : '2.0',
        license : 'GPL',
        meson_version : '>=0.49')

dialect = ''
bounds = ''
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  dialect = ['-std=f2008']
  bounds = '-fbounds-check'
  omp = '-fopenmp'
endif
if fc.get_id() == 'intel'
  dialect = ['-stand','f08']
  bounds = ['-check','bounds']
  omp = ['-qopenmp','-axAVX']
endif

add_project_arguments(dialect, language : 'fortran')
add_project_arguments(bounds,  language : 'fortran')
add_project_arguments(omp,     language : 'fortran')

cpp_std = ''
cc = meson.get_compiler('cpp')
if cc.get_id() == 'gcc'
  cpp_std = '-std=c++14'
endif
if fc.get_id() == 'intel'
  cpp_std = '-std=c++14'
endif

add_project_arguments(cpp_std, language : 'cpp')

fsrcs  = []
fapis  = []
fmain  = 'source/program_dftd.f90'
ftest  = 'source/tests_dftd.f90'
fsrcs += 'source/printout.f90'

# MCTC library
fsrcs += 'source/mctc_global.f90'
fsrcs += 'source/mctc_systools.f90'
fsrcs += 'source/mctc_readin.f90'
fsrcs += 'source/mctc_econv.f90'
fsrcs += 'source/mctc_timings.f90'
fsrcs += 'source/error.f90'
#srcs += 'source/signal.c'

# class definitions
fsrcs += 'source/class_set.f90'
fsrcs += 'source/class_molecule.f90'
fsrcs += 'source/class_param.f90'

# I/O
fsrcs += 'source/argparser.f90'
fsrcs += 'source/geometry_reader.f90'
fsrcs += 'source/disp_output.f90'

# API and interfaces
fapis += 'source/d4_calculation.f90'
fapis += 'source/d4_api.f90'

# charge model
fsrcs += 'source/eeq_model.f90'

# coordination number
fsrcs += 'source/coordination_number.f90'

# dispersion
fsrcs += 'source/dftd4.f90'
#incfil += 'include/param_ref.f'
fsrcs += 'source/dfuncpar.f90'

# C++ API example
cpp_srcs = []
cpp_srcs += 'tests/cpp_api_example.cpp'

sources = fsrcs + fapis

lapack_dep  = fc.find_library('lapack', required : true )
blas_dep    = fc.find_library('blas',   required : true )
pthread_dep = dependency('threads')

dependencies = [lapack_dep, blas_dep, pthread_dep]

incdir = include_directories('include')

# create a static library from all sources
d4lib_static = static_library(meson.project_name(), sources,
                 include_directories : incdir, pic : true)

# create a shared library as well
d4lib_shared = shared_library(meson.project_name(),
                 version : meson.project_version(),
                 dependencies : dependencies,
                 include_directories : incdir,
                 link_whole : d4lib_static,
                 link_args : '-fopenmp' )

# and link it into an executable
d4exe = executable(meson.project_name(), fmain,
                 dependencies : dependencies,
                 include_directories : incdir,
                 link_with : d4lib_static,
                 link_args : '-fopenmp' )

d4tester = executable('d4test', ftest,
                 dependencies : dependencies,
                 include_directories : incdir,
                 link_with : d4lib_static,
                 link_args : '-fopenmp' )

d4cpp = executable('dftd4_cpp', cpp_srcs,
                 dependencies : dependencies,
                 include_directories : incdir,
                 link_with : d4lib_static,
                 link_args : '-fopenmp' )

# some very basic checks to see if the executable reacts at all
test('print version',d4exe,args : '--version')
test('print help',   d4exe,args : '--help')
test('print license',d4exe,args : '--license')
test('no arguments', d4exe,should_fail : true)

# more specific tests are implemented by the tester binary
test('EEQ charges',          d4tester,args : '1')
test('Dispersion properties',d4tester,args : '2')
test('Dispersion energies',  d4tester,args : '3')
test('API tests',            d4tester,args : '4')

test('API: C++', d4cpp)
