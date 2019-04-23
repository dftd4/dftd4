DFT-D4 standalone program [![Build Status](https://travis-ci.org/awvwgk/dftd4.svg?branch=master)](https://travis-ci.org/awvwgk/dftd4)
=========================

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

Eike Caldeweyher, Christoph Bannwarth and Stefan Grimme, *J. Chem. Phys.*, **2017**, 147, 034112.
DOI: [10.1063/1.4993215](https://doi.org/10.1063/1.4993215)

Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer, Sebastian Spicher, Christoph Bannwarth and Stefan Grimme, *ChemRxiv*, **2019**, preprint. DOI: [10.26434/chemrxiv.7430216](https://doi.org/10.26434/chemrxiv.7430216.v2)

Eike Caldeweyher, Sebastian Ehlert, Andreas Hansen, Hagen Neugebauer, Sebastian Spicher, Christoph Bannwarth and Stefan Grimme,
*J. Chem. Phys.*, **2019**, 150, 154122. DOI: [https://doi.org/10.1063/1.5090222](https://doi.org/10.1063/1.5090222)


Bugs
----

please report all bugs with an example input and the used geometry,
as well as the `--verbose` output to grimme@thch.uni-bonn.de
or open an issue.

