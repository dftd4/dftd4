# 
project('dftd4', 'fortran', version : '2.0')

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

fsrcs2 = []
fsrcs2 += 'source/program_dftd.f90'
fsrcs2 += 'source/printout.f90'

# MCTC library
fsrcs2 += 'source/mctc_global.f90'
fsrcs2 += 'source/mctc_systools.f90'
fsrcs2 += 'source/mctc_readin.f90'
fsrcs2 += 'source/mctc_econv.f90'
fsrcs2 += 'source/mctc_timings.f90'
fsrcs2 += 'source/error.f90'
#srcs2 += 'source/signal.c'

# class definitions
fsrcs2 += 'source/class_set.f90'
fsrcs2 += 'source/class_molecule.f90'
fsrcs2 += 'source/class_param.f90'

# I/O
fsrcs2 += 'source/argparser.f90'
fsrcs2 += 'source/geometry_reader.f90'
fsrcs2 += 'source/disp_output.f90'

# charge model
fsrcs2 += 'source/eeq_model.f90'

# coordination number
fsrcs2 += 'source/coordination_number.f90'

# dispersion
fsrcs2 += 'source/dftd4.f90'
#incfil += 'include/param_ref.f'
fsrcs2 += 'source/dfuncpar.f90'

sources = fsrcs2

lapack_dep  = fc.find_library('lapack',  required : true)
blas_dep    = fc.find_library('blas',    required : true)
pthread_dep = dependency('threads')

dependencies = [lapack_dep, blas_dep, pthread_dep]

incdir = include_directories('include')

executable('dftd4', sources,
           dependencies : dependencies,
           include_directories : incdir,
           link_args : '-fopenmp' )
