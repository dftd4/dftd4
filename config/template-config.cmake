@PACKAGE_INIT@

set(DFTD4_WITH_OpenMP @DFTD4_WITH_OpenMP@)
set(DFTD4_WITH_API @DFTD4_WITH_API@)
set(DFTD4_WITH_API_V2 @DFTD4_WITH_API_V2@)
set(DFTD4_WITH_ILP64 @DFTD4_WITH_ILP64@)
set(DFTD4_USE_MCTCLIB @DFTD4_USE_MCTCLIB@)
set(DFTD4_USE_MULTICHARGE @DFTD4_USE_MULTICHARGE@)

enable_language("Fortran")
if(DFTD4_WITH_API)
  enable_language("C")
endif()

if(NOT TARGET "@PROJECT_NAME@::@PROJECT_NAME@")
  include("${CMAKE_CURRENT_LIST_DIR}/@PROJECT_NAME@-targets.cmake")
  include(CMakeFindDependencyMacro)

  if(NOT TARGET "OpenMP::OpenMP_Fortran" AND DFTD4_WITH_OpenMP)
    find_dependency("OpenMP")
  endif()

  if(NOT TARGET "LAPACK::LAPACK")
    find_dependency("LAPACK")
  endif()

  if(NOT TARGET "mctc-lib::mctc-lib" AND DFTD4_USE_MCTCLIB)
    find_dependency("mctc-lib")
  endif()

  if(NOT TARGET "multicharge::multicharge" AND DFTD4_USE_MULTICHARGE)
    find_dependency("multicharge")
  endif()
endif()
