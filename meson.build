# 
project('dftd4', 'fortran', 'cpp', 'c',
        version : '2.0',
        license : 'GPL',
        meson_version : '>=0.49')

doubles = ''
dialect = ''
bounds = ''
fc = meson.get_compiler('fortran')
if fc.get_id() == 'gcc'
  #doubles = '-fdefault-real-8'
  dialect = ['-std=f2008']
  bounds = '-fbounds-check'
  omp = '-fopenmp'
endif
if fc.get_id() == 'intel'
  #doubles = '-r8'
  dialect = ['-stand','f08']
  bounds = ['-check','bounds']
  omp = ['-qopenmp','-axAVX']
endif

#add_project_arguments(doubles, language : 'fortran')
add_project_arguments(dialect, language : 'fortran')
add_project_arguments(bounds,  language : 'fortran')
add_project_arguments(omp,     language : 'fortran')

cpp_std = ''
cxx = meson.get_compiler('cpp')
if cxx.get_id() == 'gcc'
  cpp_std = '-std=c++11'
endif
if cxx.get_id() == 'intel'
  cpp_std = '-std=c++11'
endif
add_project_arguments(cpp_std, language : 'cpp')

c_std = ''
cc = meson.get_compiler('c')
if cc.get_id() == 'gcc'
  c_std = '-std=c11'
endif
if cc.get_id() == 'intel'
  c_std = '-std=c11'
endif
add_project_arguments(c_std, language : 'c')


fsrcs  = []
ftest  = []
fmain  = 'source/program_dftd.f90'
ftest += 'tests/tests_dftd.f90'
ftest += 'tests/assertion.f90'
ftest += 'tests/class_molecule.f90'
ftest += 'tests/geometry_reader.f90'
ftest += 'tests/eeq_model.f90'
ftest += 'tests/dftd4.f90'
ftest += 'tests/pbc_tools.f90'
fsrcs += 'source/printout.f90'

# MCTC library
fsrcs += 'source/mctc_global.f90'
fsrcs += 'source/mctc_systools.f90'
fsrcs += 'source/mctc_readin.f90'
fsrcs += 'source/mctc_strings.f90'
fsrcs += 'source/mctc_constants.f90'
fsrcs += 'source/mctc_econv.f90'
fsrcs += 'source/mctc_param.f90'
fsrcs += 'source/mctc_timings.f90'
fsrcs += 'source/error.f90'
#srcs += 'source/signal.c'

# class definitions
fsrcs += 'source/class_set.f90'
fsrcs += 'source/class_molecule.f90'
fsrcs += 'source/class_param.f90'
fsrcs += 'source/class_wsc.f90'

# I/O
fsrcs += 'source/argparser.f90'
fsrcs += 'source/geometry_reader.f90'
fsrcs += 'source/disp_output.f90'

# API and interfaces
fsrcs += 'source/d4_calculation.f90'
fsrcs += 'source/d4_api.f90'

# periodic boundary conditions
fsrcs += 'source/pbc_tools.f90'
fsrcs += 'source/pbc.f90'

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

# C API example
c_srcs = []
c_srcs += 'tests/c_api_example.c'

sources = fsrcs

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
                 link_with : d4lib_shared)

d4c = executable('dftd4_c', c_srcs,
                 dependencies : dependencies,
                 include_directories : incdir,
                 link_with : d4lib_shared)

# some very basic checks to see if the executable reacts at all
test('Argparser: print version',d4exe,args : '--version')
test('Argparser: print help',   d4exe,args : '--help')
test('Argparser: print license',d4exe,args : '--license')
test('Argparser: no arguments', d4exe,should_fail : true)

# more specific tests are implemented by the tester binary
test('Molecule Class: axis',d4tester,args: ['class_molecule','axis'])
test('Molecule Class: MIC', d4tester,args: ['class_molecule','mic'])

test('Geometry Reader: coord 3D',d4tester,args: ['geometry_reader','coord_3d'])
test('Geometry Reader: coord 2D',d4tester,args: ['geometry_reader','coord_2d'])
test('Geometry Reader: coord 1D',d4tester,args: ['geometry_reader','coord_1d'])
test('Geometry Reader: coord 0D',d4tester,args: ['geometry_reader','coord_0d'])
test('Geometry Reader: Xmol  0D',d4tester,args: ['geometry_reader','xmol_0d'])
test('Geometry Reader: POSCAR',  d4tester,args: ['geometry_reader','poscar_3d'])

test('PBC tools: convert',   d4tester,args: ['pbc_tools','convert'])
test('PBC tools: cutoff',    d4tester,args: ['pbc_tools','cutoff'])

test('EEQ model: water',        d4tester,args: ['eeq_model','water'])
test('EEQ model: 3D Ewald',     d4tester,args: ['eeq_model','ewald'])
test('EEQ model: cell gradient',d4tester,args: ['eeq_model','cell_gradient'])

test('Dispersion: properties',    d4tester,args: ['dftd4','properties'])
test('Dispersion: energies',      d4tester,args: ['dftd4','energies'])
test('Dispersion: API',           d4tester,args: ['dftd4','api'])
test('Dispersion: API (PBC)',     d4tester,args: ['dftd4','pbc_api'])
test('Dispersion: energies (PBC)',d4tester,args: ['dftd4','pbc_disp'])
test('Dispersion: cell gradient', d4tester,args: ['dftd4','cell_gradient'])
#test('Dispersion: PBC test',      d4tester,args: ['dftd4','periodic'])

test('API: C++', d4cpp)
test('API: C', d4c)
