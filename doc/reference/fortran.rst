Fortran API
===========

The *dftd4* library seamlessly integrates with other Fortran projects via module interfaces,

.. note::

   Generally, all quantities used in the library are stored in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.

.. toctree::
   :hidden:

   Full Reference <https://dftd4.github.io/dftd4>


Handling of geometries and structure
------------------------------------

The basic infrastructure to handle molecular and periodic structures is provided by the `modular computation tool chain library <https://github.com/grimme-lab/mctc-lib>`_.
The library provides a structure type which is used to represent all geometry related informations in *dftd4*.
A structure type can be constructed from arrays or read from a file.

The constructor is provided with the generic interface ``new`` and takes an array of atomic numbers (``integer``) or element symbols (``character(len=*)``) as well as the cartesian coordinates in Bohr.
Additionally, the molecular charge and the number of unpaired electrons can be provided the ``charge`` and ``uhf`` keyword, respectively.
To create a periodic structure the lattice parameters can be passed as 3 by 3 matrix with the ``lattice`` keyword.

An example for using the constructor is given here

.. code-block:: fortran

   subroutine example
      use mctc_env, only : wp
      use mctc_io, only : structure_type, new
      implicit none
      type(structure_type) :: mol
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)

      num = [6, 1, 1, 1, 1]
      xyz = reshape([ &
        &  0.00000000000000_wp, -0.00000000000000_wp,  0.00000000000000_wp, &
        & -1.19220800552211_wp,  1.19220800552211_wp,  1.19220800552211_wp, &
        &  1.19220800552211_wp, -1.19220800552211_wp,  1.19220800552211_wp, &
        & -1.19220800552211_wp, -1.19220800552211_wp, -1.19220800552211_wp, &
        &  1.19220800552211_wp,  1.19220800552211_wp, -1.19220800552211_wp],&
        & [3, size(num)])

      call new(mol, num, xyz, charge=0.0_wp, uhf=0)

      ! ...
   end subroutine example


To interact with common input file formats for structures the ``read_structure`` procedure is available.
The file type is inferred from the name of the file automatically or if a file type hint is provided directly from the enumerator of available file types.
The ``read_structure`` routine can also use an already opened unit, but in this case the file type hint is mandatory to select the correct format to read from.

.. code-block:: fortran

   subroutine example
      use mctc_env, only : error_type
      use mctc_io, only : structure_type, read_structure, file_type
      implicit none
      type(structure_type) :: mol
      type(error_type), allocatable :: error
      character(len=:), allocatable :: input

      input = "struc.xyz"

      call read_structure(mol, input, error, file_type%xyz)
      if (allocated(error)) then
         print '(a)', error%message
         stop 1
      end if

      ! ...
   end subroutine example


The structure type as well as the error type contain only allocatable members and can therefore be used without requiring explicit deconstruction.

Certain members of the structure type should be considered immutable, like the number of atoms (``nat``), the identifiers for unique atoms (``id``) and the boundary conditions (``periodic``).
To change those specific structure parameters the structure type and all dependent objects should be reconstructed to ensure a consistent setup.
Other properties, like the geometry (``xyz``), molecular charge (``charge``), number of unpaired electrons (``uhf``) and lattice parameters (``lattice``) can be changed without requiring to reconstruct dependent objects like calculators or restart data.


Error handling
--------------

The basic error handler is an allocatable derived type, available from ``mctc_env`` as ``error_type``, which signals an error by its allocation status.

.. code-block:: fortran

   use mctc_env, only : error_type, fatal_error
   implicit none
   type(error_type), allocatable :: error

   call always_ok(error)
   if (allocated(error)) then
      print '(a)', "Unexpected failure:", error%message
   end if

   call always_failed(error)
   if (allocated(error)) then
      print '(a)', "Error:", error%message
   end if

   contains
      subroutine always_ok(error)
         type(error_type), allocatable, intent(out) :: error
      end subroutine always_ok

      subroutine always_failed(error)
         type(error_type), allocatable, intent(out) :: error

         call fatal_error(error, "Message associated with this error")
      end subroutine always_failed
   end

An unhandled error might get dropped by the next procedure call.


Performing calculations
-----------------------

An example for performing a calculation with DFT-D4 is shown below

.. code-block:: fortran

   subroutine calc_dftd4(error, mol, method, energy, model, gradient, sigma)
      use mctc_env, only : wp, error_type, fatal_error
      use mctc_io, only : structure_type
      use dftd4, only : damping_param, get_rational_damping, get_dispersion
      use dftd4_cutoff, only : realspace_cutoff
      use dftd4_model, only : dispersion_model
      use dftd4_model_d4, only : d4_model, new_d4_model
      use dftd4_model_d4s, only : d4s_model, new_d4s_model
      use dftd4_utils, only : lowercase
      implicit none

      !> Error handling
      type(error_type), allocatable, intent(out) :: error

      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      
      !> Method name for which parameters should be used
      character(len=*), intent(in) :: method
      
      !> Dispersion energy
      real(wp), intent(out) :: energy

      !> Dispersion model (D4 or D4S)
      character(len=*), intent(in), optional :: model

      !> Dispersion gradient
      real(wp), intent(out), contiguous, optional :: gradient(:, :)
      
      !> Dispersion virial
      real(wp), intent(out), contiguous, optional :: sigma(:, :)

      class(damping_param), allocatable :: param
      class(dispersion_model), allocatable :: disp

      ! Get D4 damping parameters for given method
      call get_rational_damping(method, param, s9=1.0_wp)
      if (.not.allocated(param)) then
         call fatal_error(error, "No parameters for '"//method//"' available")
         return
      end if

      ! Initialize D4/D4S model
      call new_dispersion_model(error, d4, mol)
      if (allocated(error)) return
      
      call get_dispersion(mol, disp, param, realspace_cutoff(), energy, &
         & gradient, sigma)

   end subroutine calc_dftd4


Complete Example
----------------

A minimal program using the snippets from above could look like this:

.. code-block:: fortran

   program demo
      use, intrinsic :: iso_fortran_env, only : error_unit
      use mctc_env, only : wp, error_type, fatal_error
      use mctc_io, only : new, structure_type
      use dftd4, only : damping_param, get_rational_damping, get_dispersion
      use dftd4_cutoff, only : realspace_cutoff
      use dftd4_model, only : dispersion_model, new_dispersion_model
      use dftd4_utils, only : lowercase
      implicit none
   
      !> Error handling
      type(error_type), allocatable :: error
   
      type(structure_type) :: mol
      real(wp) :: energy
      real(wp), allocatable :: gradient(:, :), sigma(:, :)
   
      call ch4(mol)      
      allocate(gradient(3, mol%nat), sigma(3, 3))
   
      call calc_dftd4(error, mol, "pbe", energy, gradient=gradient, sigma=sigma)
   
      if (allocated(error)) then
         write(error_unit, '("[Error]:", 1x, a)') error%message
         error stop
      end if
   
      write (*,'(a,f18.12)') "D4 dispersion energy (Hartree): ", energy
   
   
   contains
   
   
   subroutine ch4(mol)

      !> Molecular structure data
      type(structure_type), intent(out) :: mol
   
      real(wp), allocatable :: xyz(:, :)
      integer, allocatable :: num(:)
   
      num = [6, 1, 1, 1, 1]
      xyz = reshape([ &
        &  0.00000000000000_wp, -0.00000000000000_wp,  0.00000000000000_wp, &
        & -1.19220800552211_wp,  1.19220800552211_wp,  1.19220800552211_wp, &
        &  1.19220800552211_wp, -1.19220800552211_wp,  1.19220800552211_wp, &
        & -1.19220800552211_wp, -1.19220800552211_wp, -1.19220800552211_wp, &
        &  1.19220800552211_wp,  1.19220800552211_wp, -1.19220800552211_wp],&
        & [3, size(num)])
   
      call new(mol, num, xyz, charge=0.0_wp, uhf=0)
   
   end subroutine ch4
   
   
   subroutine calc_dftd4(error, mol, method, energy, model, gradient, sigma)
   
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
   
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      
      !> Method name for which parameters should be used
      character(len=*), intent(in) :: method
      
      !> Dispersion energy
      real(wp), intent(out) :: energy
   
      !> Dispersion model (D4 or D4S)
      character(len=*), intent(in), optional :: model
   
      !> Dispersion gradient
      real(wp), intent(out), contiguous, optional :: gradient(:, :)
      
      !> Dispersion virial
      real(wp), intent(out), contiguous, optional :: sigma(:, :)
   
      class(damping_param), allocatable :: param
      class(dispersion_model), allocatable :: disp
   
      ! Get D4 damping parameters for given method
      call get_rational_damping(method, param, s9=1.0_wp)
      if (.not.allocated(param)) then
         call fatal_error(error, "No parameters for '"//method//"' available")
         return
      end if
   
      ! Initialize D4/D4S model
      call new_dispersion_model(error, disp, mol, model)
      if (allocated(error)) return
      
      call get_dispersion(mol, disp, param, realspace_cutoff(), energy, &
         & gradient, sigma)
   
   end subroutine calc_dftd4
   
   end program demo


The program can be compiled using the following minimal `fpm.toml`.

.. code-block:: toml

   name = "d4-demo"
   
   [dependencies]
   dftd4.git = "https://github.com/dftd4/dftd4"
   dftd4.tag = "HEAD"
