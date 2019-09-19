# 
project('dftd4', 'fortran',
        version : '2.3',
        license : 'GPL',
        meson_version : '>=0.49')

conf = configuration_data()
conf.set('version', meson.project_version())
conf.set('commit', run_command(find_program('git'),'show','-s','--format=%h').stdout().strip())
conf.set('date', run_command(find_program('date'),'-I').stdout().strip())
conf.set('author', run_command(find_program('whoami')).stdout().strip())
conf.set('origin', run_command(find_program('hostname')).stdout().strip())
configure_file(input : 'source/version.f90.in', output : 'dftd4_version.fh',
  configuration : conf)

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
  if add_languages('cpp', required: false)
    cxx = meson.get_compiler('cpp')
    if cxx.get_id() == 'gcc'
      cpp_std = '-std=c++11'
    endif
    if cxx.get_id() == 'intel'
      cpp_std = '-std=c++11'
    endif
    add_project_arguments(cpp_std, language: 'cpp')
  endif
  if add_languages('c', required: false)
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
fsrcs += 'source/mctc_environment.f90'
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
fsrcs += 'source/class_results.f90'

# I/O
fsrcs += 'source/argparser.f90'
fsrcs += 'source/geometry_reader.f90'
fsrcs += 'source/disp_output.f90'

# API and interfaces
fsrcs += 'source/d4_calculation.f90'
fsrcs += 'source/c_api.f90'
fsrcs += 'source/f_api.f90'

# periodic boundary conditions
fsrcs += 'source/pbc_tools.f90'
fsrcs += 'source/pbc.f90'

# charge model
fsrcs += 'source/eeq_model.f90'

# coordination number
fsrcs += 'source/coordination_number.f90'

# dispersion
fsrcs += 'source/dftd4.f90'
fsrcs += 'source/dfuncpar.f90'

# C++ API example
cpp_srcs = []
cpp_srcs += 'tests/cpp_api_example.cpp'

# C API example
c_srcs = []
c_srcs += 'tests/c_api_example.c'

sources = fsrcs

incdir = include_directories('include')

## ========================================== ##
## LIBRARIES
## ========================================== ##
dependencies = []
dependencies_exe = []
dependencies_sha = []

if get_option('la_backend') == 'mkl'
  libmkl = [fc.find_library('pthread')]
  libmkl += fc.find_library('m')
  libmkl += fc.find_library('dl')
  if fc.get_id() == 'intel'
    libmkl_exe = [fc.find_library('mkl_intel_lp64')]
    libmkl_exe += fc.find_library('mkl_intel_thread')
  else
    libmkl_exe = [fc.find_library('mkl_gf_lp64')]
    libmkl_exe += fc.find_library('mkl_gnu_thread')
  endif
  libmkl_exe += fc.find_library('mkl_core')
  libmkl_exe += fc.find_library('iomp5')
  libmkl_sha = [fc.find_library('mkl_rt')]
  dependencies += libmkl
  dependencies_sha += libmkl_sha
  dependencies_exe += libmkl_exe
else
  dependencies += fc.find_library('blas', required : true)
  dependencies += fc.find_library('lapack', required : true)
endif

if get_option('openmp')
  if fc.get_id() == 'intel'
    add_project_arguments('-qopenmp', language : 'fortran')
    add_project_link_arguments('-qopenmp', language : 'fortran')
  else
    add_project_arguments('-fopenmp', language : 'fortran')
    add_project_link_arguments('-fopenmp', language : 'fortran')
  endif
endif

pthread_dep = dependency('threads')
dependencies += pthread_dep

# distribute dependencies for shared object and static executable
dependencies_sha += dependencies
dependencies_exe += dependencies

## ========================================== ##
## TARGETS
## ========================================== ##

# create a static library from all sources
d4lib_static = static_library(meson.project_name(), sources,
                 include_directories : incdir, pic : true)

# create a shared library as well
d4lib_shared = shared_library(meson.project_name(),
                 version : meson.project_version(),
                 dependencies : dependencies_sha,
                 include_directories : incdir,
                 link_whole : d4lib_static,
                 install : true)

# and link it into an executable
d4exe = executable(meson.project_name(), fmain,
                 dependencies : dependencies_exe,
                 include_directories : incdir,
                 link_with : d4lib_static,
                 install : true)

## ========================================== ##
## INSTALL
## ========================================== ##
install_headers('include/dftd4.h')

asciidoc = find_program('a2x', required: false)
if asciidoc.found()
   run_command(asciidoc, 'man1/dftd4.1.txt', '--format=manpage')
   install_man('man1/dftd4.1')
endif

## ========================================== ##
## TESTSUITE
## ========================================== ##
mypy = find_program('mypy', required: false)
pytest = find_program('pytest', required: false)
pylint = find_program('pylint', required: false)

d4tester = executable('d4test', ftest,
                 dependencies : dependencies_exe,
                 include_directories : incdir,
                 link_with : d4lib_static)

py_dftd4 = meson.current_source_dir() / 'python' / 'dftd4'
pypath = environment()
pypath.prepend('LD_LIBRARY_PATH', meson.current_build_dir())
pypath.prepend('PYTHONPATH', meson.current_source_dir() / 'python')

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
test('EEQ model: num. gradient',d4tester,args: ['eeq_model','numgrad'])
test('EEQ model: num. gradient (PBC)',d4tester,args: ['eeq_model','numgrad_pbc'])

test('Dispersion: properties',    d4tester,args: ['dftd4','properties'])
test('Dispersion: energies',      d4tester,args: ['dftd4','energies'])
test('Dispersion: num. gradient', d4tester,args: ['dftd4','numgrad'])
test('Dispersion: energies (PBC)',d4tester,args: ['dftd4','pbc_disp'])
test('Dispersion: num. stress',   d4tester,args: ['dftd4','numstress'])
test('Dispersion: API',           d4tester,args: ['dftd4','api'])
test('Dispersion: API (PBC)',     d4tester,args: ['dftd4','pbc_api'])
test('Dispersion: cell gradient', d4tester,args: ['dftd4','cell_gradient'])
#test('Dispersion: PBC test',      d4tester,args: ['dftd4','periodic'])

# currently we only test the API from the C/C++ side on Linux
if build_machine.system() == 'linux'
  if add_languages('cpp', required: false)
    d4cpp = executable('dftd4_cpp', cpp_srcs,
                       dependencies : dependencies,
                       include_directories : incdir,
                       link_with : d4lib_shared)
    test('API: C++', d4cpp)
  endif
  if add_languages('c', required: false)
    d4c = executable('dftd4_c', c_srcs,
                     dependencies : dependencies,
                     include_directories : incdir,
                     link_with : d4lib_shared)
    test('API: C', d4c)
  endif
endif

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
