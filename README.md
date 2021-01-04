DFT-D4 standalone program
=========================

[![Build Status](https://travis-ci.org/dftd4/dftd4.svg?branch=master)](https://travis-ci.org/dftd4/dftd4)
[![Build Status](https://github.com/dftd4/dftd4/workflows/CI/badge.svg)](https://github.com/dftd4/dftd4/actions)
[![License](https://img.shields.io/github/license/dftd4/dftd4)](https://github.com/dftd4/dftd4/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/dftd4/dftd4)](https://github.com/dftd4/dftd4/releases/latest)
[![DOI](https://zenodo.org/badge/173139980.svg)](https://zenodo.org/badge/latestdoi/173139980)

This is a minimal standalone version of DFT-D4 providing the D4(EEQ)-ATM method.


Installing
----------

To compile this version of DFT-D4 the following programs are needed
(the number in parentheses specifies the tested versions).

* `gfortran` or `ifort` compiler
* `meson` (0.53 or newer) and `ninja` (1.7 or newer) as build system
* `asciidoctor` to build the man-page

The program is build by

    $ FC=ifort meson setup build && ninja -C build

The binary is found at build/dftd4 and is ready to use.

The man-page can be build by

    $ asciidoc --doctype manpage --format manpage man1/dftd4.1.txt

By adding the directory to the `MANPATH` variable the documentation
of DFT-D4 is accessable by `man`.

`dftd4` as been successfully build using

* `ifort` 19.0.3 with the MKL as linear algebra backend
  on Manjaro Linux 18.0
* `gfortran` 8.3.0 with the MKL (19.0.3.199) as linear algebra backend
  on Manjaro Linux 18.0
* `gfortran` 8.3.0 with LAPACK (3.8.0-2) and openBLAS (0.3.6-1)
  as linear algebra backend on Manjaro Linux 18.0

`dftd4` could not be compiled with

* `gfortran` 4 or older (missing Fortran 2003 standard)


Available Package
-----------------

Statically linked binaries can be found at the [latest release page](https://github.com/dftd4/dftd4/releases/latest).
DFT-D4 is now also available for some package managers.


### Arch User Repository (AUR)

[![AUR stable](https://img.shields.io/aur/version/dftd4)](https://aur.archlinux.org/packages/dftd4/)
[![AUR git](https://img.shields.io/aur/version/dftd4-git?label=aur-git)](https://aur.archlinux.org/packages/dftd4-git/)

Get the build file from the [AUR](https://aur.archlinux.org) and use `makepkg`
(or your favorite wrapper) to build the package for `pacman`.
`dftd4` will be build with openBLAS backend, openMP and the most recent GCC.

The package build files are available [here](assets/aur) as submodules.


### Conda-Forge

[![Conda Version](https://img.shields.io/conda/vn/conda-forge/dftd4.svg)](https://anaconda.org/conda-forge/dftd4)

Installing `dftd4` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
```

Once the `conda-forge` channel has been enabled, `dftd4` can be installed with:

```
conda install dftd4
```

It is possible to list all of the versions of `dftd4` available on your platform with:

```
conda search dftd4 --channel conda-forge
```


Usage
-----

DFT-D4 is invoked by

    $ dftd4 [options] <file>

where file is a valid xyz-file (coordinates in Ångström) or a
Turbomole coord file containing only the $coord data group with
coordinates in Bohr.

More information can be obtained from the manpage or by
invoking the program intern help page with `--help`.


Examples
--------

For a general D4 calculation to obtain C6 coefficients use

    $ dftd4 coord

To calculate a dispersion correction for a PBE0 calculation use

    $ dftd4 --func pbe0 coord

Note that, the `--func` option will try to make as much sense as
possible from your input, by converting it internally to lowercase
and ignoring most dashes or year number. The input `b-p` and `BP86/def-TZVP`
are therefore equivalent.

To calculate the derivative of the dispersion energy use

    $ dftd4 --func pbe0 --grad coord

This will write a Turbomole style `gradient` file or will try to
augment an already present `gradient` file with the dispersion gradient.

For the D4(EEQ)-MBD method use

    $ dftd4 --func pbe0 --mbd coord


Citation
--------

Always cite:

Eike Caldeweyher, Christoph Bannwarth and Stefan Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
DOI: [10.1063/1.4993215](https://doi.org/10.1063/1.4993215)

Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer, Sebastian Spicher, Christoph Bannwarth and Stefan Grimme, *J. Chem Phys*, **2019**, 150, 154122.
DOI: [10.1063/1.5090222](https://doi.org/10.1063/1.5090222)
chemrxiv: [10.26434/chemrxiv.7430216](https://doi.org/10.26434/chemrxiv.7430216.v2)

Eike Caldeweyher, Jan-Michael Mewes, Sebastian Ehlert and Stefan Grimme, *Phys. Chem. Chem. Phys.*, **2020**, 22, 8499-8512.
DOI: [10.1039/D0CP00502A](https://doi.org/10.1039/D0CP00502A)
chemrxiv: [10.26434/chemrxiv.10299428](https://doi.org/10.26434/chemrxiv.10299428.v1)


License
-------

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
