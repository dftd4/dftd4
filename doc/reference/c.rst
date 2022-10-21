C API
=====

The C API bindings are provided by using the ``iso_c_binding`` intrinsic module.
Generally, objects are exported as opaque pointers and can only be manipulated within the library.
The API user is required delete all objects created in the library by using the provided deconstructor functions to avoid mamory leaks.

Overall four classes of objects are provided by the library

- error handlers (``dftd4_error``),
  used to communicate exceptional conditions and errors from the library to the user
- structure containers (``dftd4_structure``),
  used to represent the system specific information and geometry data,
  only the latter are mutable for the user
- dispersion model objects (``dftd4_model``),
  general model for calculating dispersion releated properties
- damping function objects (``dftd4_param``)
  polymorphic objects to represent the actual method parametrisation

.. note::

   Generally, all quantities provided to the library are assumed to be in `atomic units <https://en.wikipedia.org/wiki/Hartree_atomic_units>`_.


Error handling
--------------

.. c:type:: struct _dftd4_error* dftd4_error;

   Error handle class

The library provides a light error handle type (``dftd4_error``) for storing error information
The error handle requires only small overhead to construct and can only contain a single error.

The handler is represented by an opaque pointer and can only be manipulated by call from the library.
The user of those objects is required to delete the handlers again using the library provided deconstructors to avoid memory leaks.

.. c:function:: dftd4_error dftd4_new_error();

   :returns: New allocation for error handle

   Create new error handle object

.. c:function:: int dftd4_check_error(dftd4_error error);

   :param error: Error handle
   :returns: Current status of error handle, non-zero in case of error

   Check error handle status

.. c:function:: void dftd4_get_error(dftd4_error error, char* buffer, const int* buffersize);

   :param error: Error handle
   :param buffer: Allocation to store error message in
   :param buffersize: Maximum length of the buffer (optional)

   Get error message from error handle

.. c:function:: void dftd4_delete_error(dftd4_error* error);

   :param error: Error handle

   Delete error handle object


Structure data
--------------

.. c:type:: struct _dftd4_structure* dftd4_structure;

   Molecular structure data class

The structure data is used to represent the system of interest in the library.
It contains immutable system specific information like the number of atoms, the unique atom groups and the boundary conditions as well as mutable geometry data like cartesian coordinates and lattice parameters.

.. c:function:: dftd4_structure dftd4_new_structure(dftd4_error error, const int natoms, const int* numbers, const double* positions, const double* charge, const double* lattice, const bool* periodic);

   :param natoms: Number of atoms in the system
   :param numbers: Atomic numbers of all atoms [natoms]
   :param positions: Cartesian coordinates in Bohr [natoms, 3]
   :param charge: Total molecular charge (optional)
   :param lattice: Lattice parameters in Bohr [3, 3] (optional)
   :param periodic: Periodic dimension of the system [3] (optional)
   :returns: New molecular structure data handle

   Create new molecular structure data (quantities in Bohr)

.. c:function:: void dftd4_delete_structure(dftd4_structure* mol);

   :param mol: Molecular structure data handle

   Delete molecular structure data

.. c:function:: void dftd4_update_structure(dftd4_error error, dftd4_structure mol, const double* positions, const double* lattice);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param positions: Cartesian coordinates in Bohr [natoms, 3]
   :param lattice: Lattice parameters in Bohr [3, 3] (optional)

   Update coordinates and lattice parameters (quantities in Bohr)


Dispersion model
----------------

.. c:type:: struct _dftd4_model* dftd4_model;

   Dispersion model class

Instantiated for a given molecular structure type, it carries no information on the geometry but relies on the atomic species of the structure object.
Recreating a structure object requires to recreate the dispersion model as well.

.. c:function:: dftd4_model dftd4_new_d4_model(dftd4_error error, dftd4_structure mol);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :returns: New dispersion model handle

   Create new D4 dispersion model

.. c:function:: dftd4_model dftd4_custom_d4_model(dftd4_error error, dftd4_structure mol, double ga, double gc, double wf);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param ga: Charge scaling height
   :param gc: Charge scaling steepness
   :param wf: Weighting factor for coordination number interpolation
   :returns: New dispersion model handle

   Create new D4 dispersion model with custom parameters

.. c:function:: void dftd4_delete_model(dftd4_model* disp);

   :param disp: Dispersion model handle

   Delete dispersion model


Damping parameters
------------------

.. c:type:: struct _dftd4_param* dftd4_param;

   Damping parameter class

The damping parameter object determining the short-range behaviour of the dispersion correction.
Standard damping parameters like the rational damping are independent of the molecular structure and can easily be reused for several structures or easily exchanged.

.. c:function:: dftd4_param dftd4_new_rational_damping(dftd4_error error, double s6, double s8, double s9, double a1, double a2, double alp);

   :param error: Error handle
   :param s6: Scaling factor for C6 contribution
   :param s8: Scaling factor for C8 contribution
   :param s9: Scaling factor for C9 contribution
   :param a1: Scaling factor for critical radii
   :param a2: Offset distance in Bohr for critical radii
   :returns: New damping function parameter handle

   Create new rational damping parameters

.. c:function:: dftd4_param dftd4_load_rational_damping(dftd4_error error, char* method, bool mdb);

   :param error: Error handle
   :param method: Name of the method to load parameters for
   :param mbd: Use three-body specific parametrization
   :returns: New damping function parameter handle

   Load rational damping parameters from internal storage

.. c:function:: void dftd4_delete_param(dftd4_param* param);

   :param param: Damping function parameter handle

   Delete damping parameters


Calculation entrypoints
-----------------------

To evaluate dispersion energies or related properties the `dftd4_get_dispersion` procedure and similar can be used.

.. c:function:: void dftd4_get_properties(dftd4_error error, dftd4_structure mol, dftd4_model disp, double* cn, double* charges, double* c6, double* alpha);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param cn: Coordination number for all atoms [natoms]
   :param charges: Partial charges for all atoms [natoms]
   :param c6: C6 coefficients for all atom pairs [natoms, natoms]
   :param alpha: Static polarizibilities for all atoms [natoms]

   Evaluate properties related to the dispersion model

.. c:function:: void dftd4_get_dispersion(dftd4_error error, dftd4_structure mol, dftd4_model disp, dftd4_param param, double* energy, double* gradient, double* sigma);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param param: Damping function parameter handle
   :param energy: Dispersion energy
   :param gradient: Dispersion gradient [natoms, 3] (optional)
   :param sigma: Dispersion strain derivatives [3, 3] (optional)

   Evaluate the dispersion energy and its derivatives

.. c:function:: void dftd4_get_pairwise_dispersion(dftd4_error error, dftd4_structure mol, dftd4_model disp, dftd4_param param, double* pair_energy2, double* pair_energy3);

   :param error: Error handle
   :param mol: Molecular structure data handle
   :param disp: Dispersion model handle
   :param param: Damping function parameter handle
   :param energy: Pairwise additive dispersion energies
   :param energy: Pairwise non-addititive dispersion energies

   Evaluate the pairwise representation of the dispersion energy
