DFT-D4 project
==============

[![License](https://img.shields.io/github/license/dftd4/dftd4)](https://github.com/dftd4/dftd4/blob/master/COPYING.LESSER)
[![Latest Version](https://img.shields.io/github/v/release/dftd4/dftd4)](https://github.com/dftd4/dftd4/releases/latest)
[![Build Status](https://github.com/dftd4/dftd4/workflows/CI/badge.svg)](https://github.com/dftd4/dftd4/actions)
[![Documentation Status](https://readthedocs.org/projects/dftd4/badge/?version=latest)](https://dftd4.readthedocs.io/en/latest/?badge=latest)
[![docs](https://github.com/dftd4/dftd4/workflows/docs/badge.svg)](https://dftd4.github.io/dftd4/)
[![codecov](https://codecov.io/gh/dftd4/dftd4/branch/main/graph/badge.svg?token=IFtEF1Hwqj)](https://codecov.io/gh/dftd4/dftd4)

Generally Applicable Atomic-Charge Dependent London Dispersion Correction.


## Installing

A statically linked binary distribution for Linux platforms is available at the [latest release](https://github.com/dftd4/dftd4/releases/latest) tag.
Bleeding edge releases of the latest source from this repository are available on the [continuous release tag](https://github.com/dftd4/dftd4/releases/tag/bleed).


### Conda package

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/dftd4.svg?label=dftd4)](https://anaconda.org/conda-forge/dftd4)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/dftd4-python.svg?label=dftd4-python)](https://anaconda.org/conda-forge/dftd4-python)

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the [miniforge](https://github.com/conda-forge/miniforge/releases) installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

```
conda config --add channels conda-forge
```

Once the *conda-forge* channel has been enabled, this project can be installed with:

```
conda install dftd4
```

If you want to enable the Python API as well install

```
conda install dftd4-python
```

It is possible to list all of the versions available on your platform with:

```
conda search dftd4 --channel conda-forge
```

Now you are ready to use `dftd4`.


### Spack package

[![Spack](https://img.shields.io/spack/v/dftd4)](https://spack.readthedocs.io/en/latest/package_list.html#dftd4)

This project is available with [spack](https://spack.io) in its develop version.
You can install and load `dftd4` with

```
spack install dftd4
spack load dftd4
```

The Python API can be enabled by adding `+python` to the command.


### FreeBSD Port

[![FreeBSD port](https://repology.org/badge/version-for-repo/freebsd/dftd4.svg)](https://www.freshports.org/science/dftd4/)

A port for FreeBSD is available and can be installed using

```
pkg install science/dftd4
```

In case no package is available build the port using

```
cd /usr/ports/science/dftd4
make install clean
```

For more information see the [dftd4 port details](https://www.freshports.org/science/dftd4/).



### Building from Source

To build this project from the source code in this repository you need to have
a Fortran compiler supporting Fortran 2008 and one of the supported build systems:
- [meson](https://mesonbuild.com) version 0.55 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
- [cmake](https://cmake.org) version 3.14 or newer, with
  a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.10 or newer
- [fpm](https://github.com/fortran-lang/fpm) version 0.2.0 or newer

Currently this project supports GCC and Intel compilers.


#### Building with meson

To compile this version of DFT-D4 with meson the following programs are needed
(the number in parentheses specifies the tested versions).

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.55 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer
- a LAPACK / BLAS provider, like MKL or OpenBLAS

Optional dependencies are
- asciidoctor to build the manual page
- FORD to build the developer documentation
- C compiler to test the C-API and compile the Python extension module
- Python 3.6 or newer with the CFFI package installed to build the Python API

Setup a build with

```sh
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable.
To compile and run the projects testsuite use

```sh
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```sh
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.


#### Building with CMake

Alternatively, this project can be build with CMake (in this case ninja 1.10 or newer is required):

```
cmake -B _build -G Ninja -DCMAKE_INSTALL_PREFIX=$HOME/.local
```

To compile the project with CMake run

```
cmake --build _build
```

You can run the project testsuite with

```
ctest --test-dir _build --parallel --output-on-failure
```

Finally, you can install the project to the selected prefix

```
cmake --install _build
```

Note that the CMake build does not support to build the Python extension module as part of the main build.


#### Building with fpm

This project support the Fortran package manager (fpm).
Invoke fpm in the project root with

```
fpm build
```

To run the testsuite use

```
fpm test
```

You can access the ``dftd4`` program using the run subcommand

```
fpm run -- --help
```

To use ``dftd4`` for testing include it as dependency in your package manifest

```toml
[dependencies]
dftd4.git = "https://github.com/dftd4/dftd4"
```

Note that the fpm build does not support exporting the C-API, it only provides access to the standalone binary.


## Usage

DFT-D4 calculations can be performed with the `dftd4` executable.
To calculate the dispersion correction for PBE0-D4 run:

```sh
dftd4 --func pbe0 coord
```

In case you want to access the DFT-D4 results from other programs, dump the results to JSON with
(the `--noedisp` flag prevents the `.EDISP` file generation):

```sh
dftd4 --func pbe0 --json --grad --noedisp struct.xyz
```

Dispersion related properties can be calculated as well:

```sh
dftd4 --property geo.gen
```

To evaluate pairwise resolved dispersion energies use

```sh
dftd4 --pair-resolved mol.xyz
```

For an overview over all command line arguments use the `--help` argument or checkout the [`dftd4(1)`](https://github.com/dftd4/dftd4/blob/main/man/dftd4.1.adoc) manpage.


## Parameters

DFT-D4 is parametrized for plenty of density functionals.
The available parameters are listed in the [parameters.toml file](https://github.com/dftd4/dftd4/blob/main/assets/parameters.toml).

You can add new functionals using to the TOML file by adding a new subtable

```toml
[parameter.name]
reference.doi = ["<functional reference>"]
d4.bj-eeq-atm = { s8=1.0, a1=0.4, a2=5.0, doi="<parameter reference>" }
```

Those parameters are currently only used as reference and not yet usable in the library or executable.


## API access

The DFT-D4 project provides first class API support Fortran, C and Python.
Other programming languages should try to interface with to DFT-D4 via one of those three APIs.
To provide first class API support for a new language the interface specification should be available as meson build files.

The `dftd4` binary provides with the `--json` option access to all quantities available from the APIs as well.


### Fortran API

The recommended way to access the Fortran module API is by using `dftd4` as a meson subproject.
Alternatively, the project is accessible by the Fortran package manager ([fpm](https://github.com/fortran-lang/fpm)).

The complete API is available from `dftd4` module, the individual modules are available to the user as well but are not part of the public API and therefore not guaranteed to remain stable.
ABI compatibility is only guaranteed for the same minor version.

The communication with the Fortran API uses the `error_type` and `structure_type` of the modular computation tool chain library (mctc-lib) to handle errors and represent geometries, respectively.


#### Building Vasp with support for D4

To use ``dftd4`` in Vasp the compatibility layer for the 2.5.x API has to be enable with ``-Dapi_v2=true`` (meson) or ``-DWITH_API_V2=ON`` (CMake).
It is important to build ``dftd4`` with the same Fortran compiler you build Vasp with.

After you completed the installation of ``dftd4``, make sure it is findable by ``pkg-config``, you can check by running:

```
pkg-config --modversion dftd4
```

If your ``dftd4`` installation is not findable, you have to update your environment variables.
One option is to provide a module file for your ``dftd4`` installation.
The example module file below can be placed in your ``MODULEPATH`` to provide access to an installation in ``~/opt/dftd4/3.5.0``.
Retry the above comment after loading the ``dftd4`` module and adjust the module file until ``pkg-config`` finds your installation.

```lua
-- dftd4/3.5.0.lua
local name = "dftd4"
local version = "3.5.0"
local prefix = pathJoin(os.getenv("HOME"), "opt", name, version)
local libdir = "lib"  -- or lib64

whatis("Name        : " .. name)
whatis("Version     : " .. version)
whatis("Description : Generally applicable charge dependent London dispersion correction")
whatis("URL         : https://github.com/dftd4/dftd4")

prepend_path("PATH", pathJoin(prefix, "bin"))
prepend_path("MANPATH", pathJoin(prefix, "share", "man"))
prepend_path("CPATH", pathJoin(prefix, "include"))
prepend_path("LIBRARY_PATH", pathJoin(prefix, libdir))
prepend_path("LD_LIBRARY_PATH", pathJoin(prefix, libdir))
prepend_path("PKG_CONFIG_PATH", pathJoin(prefix, libdir, "pkgconfig"))
```

To enable support for D4 in Vasp add the following lines to the Makefile:

```make
CPP_OPTIONS += -DDFTD4
LLIBS       += $(shell pkg-config --libs dftd4)
INCS        += $(shell pkg-config --cflags dftd4)
```


### C API

The C API provides access to the basic Fortran objects and their most important methods to interact with them.
All Fortran objects are available as opaque `void*` in C and can only be manipulated with the correct API calls.
To evaluate a dispersion correction in C four objects are available:

1. the error handler:

   Simple error handler to carry runtime exceptions created by the library.
   Exceptions can be handled and/or transfered to the downstream error handling system by this means.

2. the molecular structure data:

   Provides a representation of the molecular structure with immutable number of atoms, atomic species, total charge and boundary conditions.
   The object provides a way to update coordinates and lattice parameters, to update immutable quantities the object has to be recreated.

3. the dispersion model:

   Instantiated for a given molecular structure type, it carries no information on the geometry but relies on the atomic species of the structure object.
   Recreating a structure object requires to recreate the dispersion model as well.

4. the damping parameters:

   Damping parameter object determining the short-range behaviour of the dispersion correction.
   Standard damping parameters like the rational damping are independent of the molecular structure and can easily be reused for several structures or easily exchanged.

The user is responsible for creating and deleting the objects to avoid memory leaks.
For convenience the type-generic macro ``dftd4_delete`` is available to free any memory allocation made in the library.


### Python API

[![LGTM](https://img.shields.io/lgtm/grade/python/g/dftd4/dftd4.svg)](https://lgtm.com/projects/g/dftd4/dftd4/context:python)

The Python API is disabled by default and can be built in-tree or out-of-tree.
The in-tree build is mainly meant for end users and packages.
To build the Python API with the normal project set the `python` option in the configuration step with

```sh
meson setup _build -Dpython=true -Dpython_version=$(which python3)
```

The Python version can be used to select a different Python version, it defaults to `'python3'`.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Proceed with the build as described before and install the projects to make the Python API available in the selected prefix.

For the out-of-tree build see the instructions in the [`python`](https://github.com/dftd4/dftd4/blob/main/python) directory.


## Citation

Always cite:

Eike Caldeweyher, Christoph Bannwarth and Stefan Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
DOI: [10.1063/1.4993215](https://doi.org/10.1063/1.4993215)

Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer, Sebastian Spicher, Christoph Bannwarth and Stefan Grimme, *J. Chem Phys*, **2019**, 150, 154122.
DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)
chemrxiv: [10.26434/chemrxiv.7430216](https://doi.org/10.26434/chemrxiv.7430216.v2)

Eike Caldeweyher, Jan-Michael Mewes, Sebastian Ehlert and Stefan Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 8499-8512.
DOI: [10.1039/D0CP00502A](https://doi.org/10.1039/D0CP00502A)
chemrxiv: [10.26434/chemrxiv.10299428](https://doi.org/10.26434/chemrxiv.10299428.v1)


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the Lesser GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
Lesser GNU General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
Lesser GNU General Public license, shall be licensed as above, without any
additional terms or conditions.
