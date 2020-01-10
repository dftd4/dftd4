DFT-D4 standalone program
=========================

[![Build Status](https://travis-ci.org/dftd4/dftd4.svg?branch=master)](https://travis-ci.org/dftd4/dftd4)
[![License](https://img.shields.io/github/license/dftd4/dftd4)](https://github.com/dftd4/dftd4/blob/master/COPYING)
[![Latest Version](https://img.shields.io/github/v/release/dftd4/dftd4)](https://github.com/dftd4/dftd4/releases/latest)

Copied from
https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dftd4

This is a minimal standalone version of DFT-D4 providing the
D4(EEQ)-ATM and D4(EEQ)-MBD methods.

Installing
----------

To compile this version of DFT-D4 the following programs are needed
(the number in parentheses specifies the tested versions).

* `gfortran` (8.2.1) or `ifort` (17.0.7 or 18.0.3) compiler
* `meson` (0.49.0) and `ninja` (1.8.2) as build system
* `asciidoc` (8.6.10) to build the man-page

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

DFT-D4 is now also available for some package managers.

### Arch User Repository (AUR)

[![AUR stable](https://img.shields.io/aur/version/dftd4)](https://aur.archlinux.org/packages/dftd4/)
[![AUR git](https://img.shields.io/aur/version/dftd4-git?label=aur-git)](https://aur.archlinux.org/packages/dftd4-git/)

Clone the AUR repository and use `makepkg` (or your favorite wrapper) to build
the package for `pacman`.
You can choose between the current stable and the master version from this
repository.
The build will use openBLAS and the most up-to-date GCC.

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

Eike Caldeweyher, Jan-Michael Mewes, Sebastian Ehlert and Stefan Grimme, *ChemRxiv*, **2019**, preprint.
chemrxiv: [10.26434/chemrxiv.10299428](https://doi.org/10.26434/chemrxiv.10299428.v1)

Bugs
----

please report all bugs with an example input and the used geometry,
as well as the `--verbose` output to [Stefan Grimme](mailto:grimme@thch.uni-bonn.de)
or open an issue.
