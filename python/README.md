# DFT-D4 Python wrapper

This is the Python-side of the wrapper around the `dftd4` library.
It provides access to all features of the binary and additional ones
not accessible from the binary. Also there are some functions provided
to manipulate data produced in `dftd4`.

## Usage

For the everyday use we recommend to use `dftd4` together with the
Atomic Simulation Model (ASE) like

```python
from ase.io import read
from dftd4 import D4_model
mol = read('coord')
mol.set_calculator(D4_model(xc='pbe0'))
energy = mol.get_potential_energy()
forces = mol.get_forces()
```

To use `dftd4` as a dispersion correction rather than a standalone,
it can be chained with another calculator.

```python
from ase.io import read
from gpaw import GPAW, PW
from dftd4 import D4_model
mol = read('POSCAR')
calc = D4_model(xc='PBE', calc=GPAW(xc='PBE', mode=PW(300), txt='gpaw.out'))
mol.set_calculator(calc)
energy = mol.get_potential_energy()
forces = mol.get_forces()
```

Additionally the dispersion calculators provide the possibility to
return properties like polarizabilities or C6-coefficients.
For example one could use the utilities submodule to extrapolate higher
order dispersion coefficients from the calculated C6-coefficients.

```python
from dftd4 import D4_model
from dftd4.utils import extrapolate_c_n_coeff
from numpy import zeros
# get atoms from somewhere
calc = D4_model(energy=False, forces=False)
c6 = calc.get_property('c6_coefficients', atoms=atoms)
c8 = zeros((len(atoms), len(atoms)))
numbers = atoms.get_atomic_numbers()
for i in range(len(atoms)):
    for j in range(len(atoms)):
        c8[i, j] = extrapolate_c_n_coeff(c6[i, j], numbers[i], numbers[j], 8)
```

## Technical Details

We provide a multilayered API for `dftd4` starting from the Fortran side
with a set of standalone calculator subroutines in `source/d4_calculation.f90`,
these calculators are exposed by using the `iso_c_binding` module as C-API
(see `source/c_api.f90` and `include/dftd4.h`).
All non-Fortran compatible languages (that makes all languages except for Fortran),
should use this C-API to interface with `dftd4`.

To use the C-API in Python we decided to define the interface on the Python
side rather than on the Fortran side using the `ctypes` module.
This interface (or Python-API) hides some of the implementation dependent details
from the user and provides access to all main functions implemented in `dftd4`.

Using the Python-API we wrapped everything up for the Atomic Simulation Environment
(ASE) by implementing a `Calculator` that handles the conversion between the
Atoms objects and the bundle of ndarrays required for Python-API.

For Python users we recommend using the ASE Calculators which give access to a
feature-rich environment for computational chemistry.
For more information on ASE we refer to their [detailed documentation](https://wiki.fysik.dtu.dk/ase/).

If you want to interface your Python program with `dftd4` and are afraid to add
a such heavy module as ASE to your dependencies, you can interface directly
to the Python-API which lives in `dftd4.interface` and does only require
`numpy` and `ctypes`.

Coming from every other C-compatible language we would recommend to wrap
the C-API in a similar way like we did for Python.
    
## Citation

Eike Caldeweyher, Christoph Bannwarth and Stefan Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
DOI: [10.1063/1.4993215](https://doi.org/10.1063/1.4993215)

Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer, Sebastian Spicher, Christoph Bannwarth and Stefan Grimme, *J. Chem Phys*, **2019**, 150, 154122. DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)
