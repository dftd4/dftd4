Using DFT-D4 in Vasp
====================

To use the D4 dispersion correction in Vasp you need to a compile the ``dftd4`` package using the same Fortran compiler as used for Vasp and enable the API compatibility needed for linking the Vasp-side interface for ``dftd4``.
Checkout the instructions for installing ``dftd4`` at :ref:`installation`.
To use ``dftd4`` in Vasp the compatibility layer for the 2.5.x API has to be enable with ``-Dapi_v2=true`` (meson) or ``-DWITH_API_V2=ON`` (CMake).

.. important::

   It is important to build ``dftd4`` with the same Fortran compiler you build Vasp with.

After you completed the installation of ``dftd4``, make sure it is findable by ``pkg-config``, you can check by running:

.. code-block:: text

   pkg-config --modversion dftd4


If your ``dftd4`` installation is not findable, you have to update your environment variables.
One option is to provide a module file for your ``dftd4`` installation.
The example module file below can be placed in your ``MODULEPATH`` to provide access to an installation in ``~/opt/dftd4/3.6.0``.
Retry the above comment after loading the ``dftd4`` module and adjust the module file until ``pkg-config`` finds your installation.

.. code-block:: lua

   -- dftd4/3.6.0.lua
   local name = "dftd4"
   local version = "3.6.0"
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


To enable support for D4 in Vasp add the following lines to the Makefile:

.. code-block:: make

   CPP_OPTIONS += -DDFTD4
   LLIBS       += $(shell pkg-config --libs dftd4)
   INCS        += $(shell pkg-config --cflags dftd4)
