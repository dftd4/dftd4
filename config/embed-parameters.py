#!/usr/bin/env python3
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
# You should have received a copy of the GNU Lesser General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

"""Generate a Fortran module containing an embedded TOML document."""

from pathlib import Path
import sys


def quote(value: str) -> str:
    """Quote one TOML fragment as a Fortran character literal."""
    return '"' + value.replace('"', '""') + '"'


def write_module(source: Path, target: Path) -> None:
    lines = source.read_text(encoding="utf-8").splitlines()
    max_line_length = max((len(line) for line in lines), default=0)
    with target.open("w", encoding="utf-8", newline="\n") as output:
        output.write("! This file is part of dftd4.\n")
        output.write("! SPDX-Identifier: LGPL-3.0-or-later\n")
        output.write("!\n")
        output.write(
            "! dftd4 is free software: you can redistribute it and/or modify it under\n"
        )
        output.write(
            "! the terms of the Lesser GNU General Public License as published by\n"
        )
        output.write(
            "! the Free Software Foundation, either version 3 of the License, or\n"
        )
        output.write("! (at your option) any later version.\n")
        output.write("!\n")
        output.write("! dftd4 is distributed in the hope that it will be useful,\n")
        output.write(
            "! but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        )
        output.write(
            "! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
        )
        output.write("! Lesser GNU General Public License for more details.\n")
        output.write("!\n")
        output.write(
            "! You should have received a copy of the GNU Lesser General Public License\n"
        )
        output.write(
            "! along with dftd4.  If not, see <https://www.gnu.org/licenses/>.\n"
        )
        output.write("!\n")
        output.write("! Generated from assets/parameters.toml; do not edit directly.\n")
        output.write("module dftd4_parameters\n")
        output.write("   implicit none\n")
        output.write("   private\n")
        output.write("\n")
        output.write("   public :: get_embedded_parameters\n")
        output.write("\n")
        output.write(f"   integer, parameter :: nlines = {len(lines)}\n")
        output.write(f"   integer, parameter :: max_line_length = {max_line_length}\n")
        output.write(
            "   character(len=max_line_length), parameter :: embedded_parameters(nlines) = [ &\n"
        )
        output.write("      & character(len=max_line_length) :: &\n")

        for index, line in enumerate(lines):
            fragments = [line[pos : pos + 96] for pos in range(0, len(line), 96)] or [
                ""
            ]
            for fragment_index, fragment in enumerate(fragments):
                last_fragment = fragment_index == len(fragments) - 1
                output.write(f"      & {quote(fragment)}")
                if not last_fragment:
                    output.write("// &\n")
                elif index < len(lines) - 1:
                    output.write(", &\n")
                else:
                    output.write(" &\n")

        output.write("      & ]\n\n")
        output.write("contains\n\n")
        output.write("pure function get_embedded_parameters() result(string)\n")
        output.write("   character(len=:), allocatable :: string\n")
        output.write("   integer :: i, length, position, line_length\n\n")
        output.write("   length = sum(len_trim(embedded_parameters)) + nlines - 1\n")
        output.write("   allocate(character(len=length) :: string)\n")
        output.write("   position = 1\n")
        output.write("   do i = 1, nlines\n")
        output.write("      line_length = len_trim(embedded_parameters(i))\n")
        output.write("      if (line_length > 0) then\n")
        output.write(
            "         string(position:position + line_length - 1) = embedded_parameters(i)(:line_length)\n"
        )
        output.write("         position = position + line_length\n")
        output.write("      end if\n")
        output.write("      if (i < nlines) then\n")
        output.write('         string(position:position) = new_line("a")\n')
        output.write("         position = position + 1\n")
        output.write("      end if\n")
        output.write("   end do\n")
        output.write("end function get_embedded_parameters\n\n")
        output.write("end module dftd4_parameters\n")


if len(sys.argv) != 3:
    raise SystemExit(f"usage: {sys.argv[0]} INPUT OUTPUT")

write_module(Path(sys.argv[1]), Path(sys.argv[2]))
