DFT-D4 Python API
-----------------

Python interface for the generally applicable atomic-charge dependent London dispersion correction, DFT-D4.
This Python project is targeted at developers who want to interface their project via Python with ``dftd4``.

This interface provides access to the C-API of ``dftd4`` via the CFFI module.
The low-level CFFI interface is available in the ``dftd4.libdftd4`` module and only required for implementing other interfaces.
A more pythonic interface is provided in the ``dftd4.interface`` module which can be used to build more specific interfaces.

.. code:: python

   >>> from dftd4.interface import DampingParam, DispersionModel
   >>> import numpy as np
   >>> numbers = np.array([1, 1, 6, 5, 1, 15, 8, 17, 13, 15, 5, 1, 9, 15, 1, 15])
   >>> positions = np.array([  # Coordinates in Bohr
   ...     [+2.79274810283778, +3.82998228828316, -2.79287054959216],
   ...     [-1.43447454186833, +0.43418729987882, +5.53854345129809],
   ...     [-3.26268343665218, -2.50644032426151, -1.56631149351046],
   ...     [+2.14548759959147, -0.88798018953965, -2.24592534506187],
   ...     [-4.30233097423181, -3.93631518670031, -0.48930754109119],
   ...     [+0.06107643564880, -3.82467931731366, -2.22333344469482],
   ...     [+0.41168550401858, +0.58105573172764, +5.56854609916143],
   ...     [+4.41363836635653, +3.92515871809283, +2.57961724984000],
   ...     [+1.33707758998700, +1.40194471661647, +1.97530004949523],
   ...     [+3.08342709834868, +1.72520024666801, -4.42666116106828],
   ...     [-3.02346932078505, +0.04438199934191, -0.27636197425010],
   ...     [+1.11508390868455, -0.97617412809198, +6.25462847718180],
   ...     [+0.61938955433011, +2.17903547389232, -6.21279842416963],
   ...     [-2.67491681346835, +3.00175899761859, +1.05038813614845],
   ...     [-4.13181080289514, -2.34226739863660, -3.44356159392859],
   ...     [+2.85007173009739, -2.64884892757600, +0.71010806424206],
   ... ])
   >>> model = DispersionModel(numbers, positions)
   >>> res = model.get_dispersion(DampingParam(method="scan"), grad=False)
   >>> res.get("energy")  # Results in atomic units
   -0.005328888532435093
   >>> res.update(**model.get_properties())  # also allows access to properties
   >>> res.get("c6 coefficients")[0, 0]
   1.5976689760849156
   >>> res.get("polarizibilities")
   array([ 1.97521745,  1.48512704,  7.33564674, 10.28920458,  1.99973802,
          22.85298573,  6.65877552, 15.39410319, 22.73119177, 22.86303028,
          14.56038118,  1.4815783 ,  3.91266859, 25.8236368 ,  1.93444627,
          23.02494331])


Additional features
~~~~~~~~~~~~~~~~~~~

The ``dftd4.parameters`` module becomes available if a TOML parser is available, either `tomlkit <https://github.com/sdispater/tomlkit>`_ or `toml <https://github.com/uiri/toml>`_ can be used here.
The returned dict can be used to supply parameters to the constructor of the ``DampingParam`` object, only the ``s6``, ``s8``, ``s9``, ``a1``, ``a2`` and ``alp`` entries will be used, the remaining entries are meta data describing the damping parameters.

.. code-block:: python

   >>> from dftd4.parameters import get_damping_param
   >>> get_damping_param("r2scan")
   {'s6': 1.0, 's9': 1.0, 'alp': 16.0, 'damping': 'bj', 'mbd': 'approx-atm', 's8': 0.6018749, 'a1': 0.51559235, 'a2': 5.77342911, 'doi': '10.1063/5.0041008'}
   >>> get_damping_param("b97m")
   {'s6': 1.0, 's9': 1.0, 'alp': 16.0, 'damping': 'bj', 'mbd': 'approx-atm', 's8': 0.6633, 'a1': 0.4288, 'a2': 3.9935, 'doi': '10.1002/jcc.26411'}


QCSchema Integration
~~~~~~~~~~~~~~~~~~~~

This Python API natively understands QCSchema and the `QCArchive infrastructure <http://docs.qcarchive.molssi.org>`_.
If the QCElemental package is installed the ``dftd4.qcschema`` module becomes importable and provides the ``run_qcschema`` function.

.. code:: python

   >>> from dftd4.qcschema import run_qcschema
   >>> import qcelemental as qcel
   >>> atomic_input = qcel.models.AtomicInput(
   ...     molecule = qcel.models.Molecule(
   ...         symbols = ["O", "H", "H"],
   ...         geometry = [
   ...             0.00000000000000,  0.00000000000000, -0.73578586109551,
   ...             1.44183152868459,  0.00000000000000,  0.36789293054775,
   ...            -1.44183152868459,  0.00000000000000,  0.36789293054775
   ...         ],
   ...     ),
   ...     driver = "energy",
   ...     model = {
   ...         "method": "TPSS-D4",
   ...     },
   ...     keywords = {},
   ... )
   ...
   >>> atomic_result = run_qcschema(atomic_input)
   >>> atomic_result.return_result
   -0.0002667885779142513


ASE Integration
~~~~~~~~~~~~~~~

To integrate with `ASE <https://wiki.fysik.dtu.dk/ase/>`_ this interface implements an ASE Calculator.
The ``DFTD4`` calculator becomes importable if an ASE installation is available.

.. code:: python

   >>> from ase.build import molecule
   >>> from dftd4.ase import DFTD4
   >>> atoms = molecule('H2O')
   >>> atoms.calc = DFTD4(method="TPSS")
   >>> atoms.get_potential_energy()
   -0.007310393443152083
   >>> atoms.calc.set(method="PBE")
   {'method': 'PBE'}
   >>> atoms.get_potential_energy()
   -0.005358475432239303
   >>> atoms.get_forces()
   array([[-0.        , -0.        ,  0.00296845],
          [-0.        ,  0.00119152, -0.00148423],
          [-0.        , -0.00119152, -0.00148423]])

To use the ``DFTD4`` calculator as dispersion correction the calculator can be combined using the `SumCalculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/mixing.html>`_ from the ``ase.calculators.mixing`` module.

.. code:: python

   >>> from ase.build import molecule
   >>> from ase.calculators.mixing import SumCalculator
   >>> from ase.calculators.nwchem import NWChem
   >>> from dftd4.ase import DFTD4
   >>> atoms = molecule('H2O')
   >>> atoms.calc = SumCalculator([DFTD4(method="PBE"), NWChem(xc="PBE")])

For convenience ``DFTD4`` allows to combine itself with another calculator by using the ``add_calculator`` method which returns a SumCalculator:

.. code:: python

   >>> from ase.build import molecule
   >>> from ase.calculators.emt import EMT
   >>> from dftd4.ase import DFTD4
   >>> atoms = molecule("C60")
   >>> atoms.calc = DFTD4(method="pbe").add_calculator(EMT())
   >>> atoms.get_potential_energy()
   6.348142387048062
   >>> [calc.get_potential_energy() for calc in atoms.calc.calcs]
   [-6.015477436263984, 12.363619823312046]

The individual contributions are available by iterating over the list of calculators in ``calc.calcs``.
Note that ``DFTD4`` will always place itself as first calculator in the list.


Installing
~~~~~~~~~~

.. image:: https://img.shields.io/conda/vn/conda-forge/dftd4-python.svg
   :alt: Conda Version
   :target: https://anaconda.org/conda-forge/dftd4-python

This project is packaged for the *conda* package manager and available on the *conda-forge* channel.
To install the *conda* package manager we recommend the `miniforge <https://github.com/conda-forge/miniforge/releases>`_ installer.
If the *conda-forge* channel is not yet enabled, add it to your channels with

.. code:: sh

   conda config --add channels conda-forge

Once the *conda-forge* channel has been enabled, this project can be installed with:

.. code:: sh

   conda install dftd4-python

Now you are ready to use ``dftd4``, check if you can import it with

.. code:: python

   >>> import dftd4
   >>> from dftd4.libdftd4 import get_api_version
   >>> get_api_version()
   '3.2.0'


Building the extension module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To perform an out-of-tree build some version of ``dftd4`` has to be available on your system and preferably findable by ``pkg-config``.
Try to find a ``dftd4`` installation you build against first with

.. code:: sh

   pkg-config --modversion dftd4

Adjust the ``PKG_CONFIG_PATH`` environment variable to include the correct directories to find the installation if necessary.


Using pip
^^^^^^^^^

This project support installation with pip as an easy way to build the Python API.

- C compiler to build the C-API and compile the extension module (the compiler name should be exported in the ``CC`` environment variable)
- Python 3.6 or newer
- The following Python packages are required additionally

  - `cffi <https://cffi.readthedocs.io/>`_
  - `numpy <https://numpy.org/>`_
  - `pkgconfig <https://pypi.org/project/pkgconfig/>`_ (setup only)

Make sure to have your C compiler set to the ``CC`` environment variable

.. code:: sh

   export CC=gcc

Install the project with pip

.. code:: sh

   pip install .

To install extra dependencies as well use

.. code:: sh

   pip install '.[parameters,ase,qcschema]'



Using meson
^^^^^^^^^^^

This directory contains a separate meson build file to allow the out-of-tree build of the CFFI extension module.
The out-of-tree build requires

- C compiler to build the C-API and compile the extension module
- `meson <https://mesonbuild.com>`_ version 0.53 or newer
- a build-system backend, *i.e.* `ninja <https://ninja-build.org>`_ version 1.7 or newer
- Python 3.6 or newer with the `CFFI <https://cffi.readthedocs.io/>`_ package installed

Setup a build with

.. code:: sh

   meson setup _build -Dpython_version=3

The Python version can be used to select a different Python version, it defaults to ``'3'``.
Python 2 is not supported with this project, the Python version key is meant to select between several local Python 3 versions.

Compile the project with

.. code:: sh

   meson compile -C _build

The extension module is now available in ``_build/dftd4/_libdftd4.*.so``.
You can install as usual with

.. code:: sh

   meson configure _build --prefix=/path/to/install
   meson install -C _build
