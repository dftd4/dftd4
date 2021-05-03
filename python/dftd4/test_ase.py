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

from dftd4.ase import DFTD4
from ase.build import molecule
from ase.calculators.emt import EMT
from pytest import approx, raises
import numpy as np


def test_ase_scand4():
    thr = 1.0e-6

    forces = np.array(
        [
            [-7.90552684e-21, -1.15811595e-19, -2.80061133e-05],
            [+4.61502216e-20, +4.26735028e-04, +4.94269127e-04],
            [+2.28682834e-19, -4.26735028e-04, +4.94269127e-04],
            [+1.39725405e-20, -9.71924142e-20, -7.95205193e-04],
            [-1.08042249e-04, -8.23929519e-05, +2.31098749e-05],
            [+1.08042249e-04, -8.23929519e-05, +2.31098749e-05],
            [+1.08042249e-04, +8.23929519e-05, +2.31098749e-05],
            [-1.08042249e-04, +8.23929519e-05, +2.31098749e-05],
            [+1.07391220e-20, -4.98420762e-05, -1.28883224e-04],
            [+5.97028977e-21, +4.98420762e-05, -1.28883224e-04],
        ]
    )

    atoms = molecule("methylenecyclopropane")
    atoms.calc = DFTD4(method="SCAN")

    assert approx(atoms.get_potential_energy(), abs=thr) == -0.021665446836610567
    assert approx(atoms.get_forces(), abs=thr) == forces

    atoms.calc = DFTD4(method="SCAN").add_calculator(EMT())
    assert approx(atoms.get_potential_energy(), abs=thr) == 3.6624398683434225
    energies = [calc.get_potential_energy() for calc in atoms.calc.calcs]
    assert approx(energies, abs=thr) == [-0.021665446836610563, 3.684105315180033]


def test_ase_tpssd4():
    thr = 1.0e-6

    forces = np.array(
        [
            [-8.27697238e-04, -1.74116189e-02, +9.50315402e-05],
            [-3.00615825e-04, +3.51410800e-04, -3.12339518e-03],
            [-9.55691674e-04, -4.05325537e-03, +7.33934018e-05],
            [+2.70525425e-04, -1.91784287e-03, -6.13796425e-04],
            [+1.56444473e-02, +5.71643192e-03, +6.29706049e-04],
            [-1.43399413e-02, +7.36397630e-03, +7.87584027e-04],
            [+4.41551907e-03, +4.04705396e-04, +1.42098826e-03],
            [-4.17670039e-03, +1.57923335e-03, +1.60488604e-03],
            [+4.66065256e-03, +7.97764912e-05, -3.60249901e-05],
            [-4.95386767e-03, +1.38115911e-03, -1.63079386e-05],
            [+5.36717422e-03, +2.78165913e-03, -4.68647341e-04],
            [-4.80380448e-03, +3.72436466e-03, -3.53417437e-04],
        ]
    )

    atoms = molecule("C2H6CHOH")
    atoms.calc = DFTD4(method="TPSS")

    assert approx(atoms.get_potential_energy(), abs=thr) == -0.24206732765720423
    assert approx(atoms.get_forces(), abs=thr) == forces

    atoms.calc = DFTD4(method="TPSS").add_calculator(EMT())
    assert approx(atoms.get_potential_energy(), abs=thr) == 4.864016486351274
    energies = [calc.get_potential_energy() for calc in atoms.calc.calcs]
    assert approx(energies, abs=thr) == [-0.24206732765720396, 5.106083814008478]
