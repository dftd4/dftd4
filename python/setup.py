# This file is part of dftd4.
#
# Copyright (C) 2019 Sebastian Ehlert
#
# dftd4 is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# dftd4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with dftd4.  If not, see <https://www.gnu.org/licenses/>.

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name="dftd4",
      version="2.4.0",
      author="Sebastian Ehlert",
      author_email="ehlert@thch.uni-bonn.de",
      description="Wrapper for the DFT-D4 program",
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords="dispersion",
      url="https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/dftd4",
      packages=find_packages(),
      install_requires=['ase'],
      classifiers=[
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: LGPL3",
          "Operating System :: Linux",
      ],
      python_requires=">=3.5",
)
