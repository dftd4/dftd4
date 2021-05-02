# This file is part of dftd4.
# SPDX-Identifier: LGPL-3.0-or-later
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# Lesser GNU General Public License for more details.
#
# You should have received a copy of the Lesser GNU General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.
"""
FFI builder module for dftd4 for usage from meson and from setup.py.

Since meson has the full knowledge about the build, it will handle
the generation of the C definitions in the meson.build file rather
than in the FFI builder. This allows to correctly keep track of
dependencies and updates in the build process.

For setup.py we have to do the preprocessing ourselves here, this
requires us to use the C compiler to preprocess the header file
of dftd4 because the CFFI C parser cannot handle certain C
preprocessor constructs. Also, we cannot rely on an external build
system fixing dependencies for us and therefore have to find those
ourselves using pkg-config.
"""

import os
import cffi

if __name__ == "__main__":
    import sys

    kwargs = dict(libraries=["dftd4"])

    header_file = sys.argv[1]
    module_name = sys.argv[2]

    with open(header_file) as f:
        cdefs = f.read()
else:
    import subprocess
    import pkgconfig

    if not pkgconfig.exists("dftd4"):
        raise Exception("Unable to find pkg-config package 'dftd4'")
    if pkgconfig.installed("dftd4", "< 3.0"):
        raise Exception(
            "Installed 'dftd4' version is too old, 3.0 or newer is required"
        )

    kwargs = pkgconfig.parse("dftd4")

    if "CC" not in os.environ:
        raise Exception("This build script requires to set a C compiler in CC")
    cc = os.environ["CC"]

    cflags = pkgconfig.cflags("dftd4").split()

    module_name = "dftd4._libdftd4"

    p = subprocess.Popen(
        [cc, *cflags, "-E", "-"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = p.communicate(b'#include "dftd4.h"')

    cdefs = out.decode()

ffibuilder = cffi.FFI()
ffibuilder.set_source(module_name, '#include "dftd4.h"', **kwargs)
ffibuilder.cdef(cdefs)

if __name__ == "__main__":
    ffibuilder.distutils_extension(".")
