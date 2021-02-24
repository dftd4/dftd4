DFT-D4 Python API
-----------------

Python interface for the generally applicable atomic-charge dependent London dispersion correction, DFT-D4.


Building the extension module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This directory contains a separate meson build file to allow the out-of-tree build of the CFFI extension module.
To perform an out-of-tree build some version of ``dftd4`` has to be available on your system and preferably findable by ``pkg-config``.
Try to find a ``dftd4`` installation you build against first with

.. code:: sh

   pkg-config --modversion dftd4

Adjust the ``PKG_CONFIG_PATH`` environment variable to include the correct directories to find the installation if necessary.
The out-of-tree build requires
- C compiler to test the C-API and compile the extension module
- `meson <https://mesonbuild.com>`_ version 0.53 or newer
- a build-system backend, *i.e.* `ninja <https://ninja-build.org>`_ version 1.7 or newer
- Python 3.6 or newer with the CFFI package installed

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

Alternatively, copy the compiled extension module into the Python module tree with

.. code:: sh

   cp _build/dftd4/_libdftd4.*.so dftd4/

Copying the extension module usually works best if the ``dftd4`` version was linked statically into the extension module.
The copying strategy is *only* recommended if you plan to develop locally at the Python API.

Finally, you can install this project with ``pip`` in development mode

.. code:: sh

   pip install -e .
