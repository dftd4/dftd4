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

import numpy as np
import pytest
from pytest import approx, raises

try:
    import pyscf
    from pyscf import lib, gto, scf
    import dftd4.pyscf as disp
except ModuleNotFoundError:
    pyscf = None


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_energy_r2scan_d4():

    mol = gto.M(
        atom="""
             C   -0.755422531  -0.796459123  -1.023590391
             C    0.634274834  -0.880017014  -1.075233285
             C    1.406955202   0.199695367  -0.653144334
             C    0.798863737   1.361204515  -0.180597909
             C   -0.593166787   1.434312023  -0.133597923
             C   -1.376239198   0.359205222  -0.553258516
             I   -1.514344238   3.173268101   0.573601106
             H    1.110906949  -1.778801728  -1.440619836
             H    1.399172302   2.197767355   0.147412751
             H    2.486417780   0.142466525  -0.689380574
             H   -2.454252250   0.422581120  -0.512807958
             H   -1.362353593  -1.630564523  -1.348743149
             S   -3.112683203   6.289227834   1.226984439
             H   -4.328789697   5.797771251   0.973373089
             C   -2.689135032   6.703163830  -0.489062886
             H   -1.684433029   7.115457372  -0.460265708
             H   -2.683867206   5.816530502  -1.115183775
             H   -3.365330613   7.451201412  -0.890098894
             """
    )

    d4 = disp.DFTD4Dispersion(mol, xc="r2SCAN")
    assert d4.kernel()[0] == approx(-0.005001101058518388, abs=1.0e-7)


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_gradient_b97m_d4():
    mol = gto.M(
        atom="""
             H    0.002144194   0.361043475   0.029799709
             C    0.015020592   0.274789738   1.107648016
             C    1.227632658   0.296655040   1.794629427
             C    1.243958826   0.183702791   3.183703934
             C    0.047958213   0.048915002   3.886484583
             C   -1.165135654   0.026954348   3.200213281
             C   -1.181832083   0.139828643   1.810376587
             H    2.155807907   0.399177037   1.249441585
             H    2.184979344   0.198598553   3.716170761
             H    0.060934662  -0.040672756   4.964014252
             H   -2.093220602  -0.078628959   3.745125056
             H   -2.122845437   0.123257119   1.277645797
             Br  -0.268325907  -3.194209024   1.994458950
             C    0.049999933  -5.089197474   1.929391171
             F    0.078949601  -5.512441335   0.671851563
             F    1.211983937  -5.383996300   2.498664481
             F   -0.909987405  -5.743747328   2.570721738
             """
    )
    ref = np.array(
        [
            [+6.02987248e-07, +1.18181692e-05, -2.11659178e-05],
            [+3.77083487e-07, +4.21255367e-05, -3.65576556e-05],
            [+3.71749233e-05, +4.38986750e-05, -1.64037320e-05],
            [+3.79004788e-05, +4.09262181e-05, +2.57427629e-05],
            [+1.49281462e-06, +3.63132380e-05, +4.66732244e-05],
            [-3.45592945e-05, +3.46256250e-05, +2.53829747e-05],
            [-3.48859913e-05, +3.74107269e-05, -1.56473785e-05],
            [+2.00543104e-05, +1.15042699e-05, -9.90469697e-06],
            [+1.99879228e-05, +9.25641402e-06, +1.21976769e-05],
            [+1.10396127e-06, +7.69249859e-06, +2.38607706e-05],
            [-1.86258815e-05, +7.79467748e-06, +1.29284817e-05],
            [-1.87883833e-05, +9.46661745e-06, -9.65731010e-06],
            [-2.38952311e-05, -1.10356928e-04, -2.28127181e-05],
            [+4.05848507e-07, -5.94239995e-05, -6.36138164e-06],
            [+2.78030538e-06, -3.80326610e-05, -1.91595254e-05],
            [+1.91553258e-05, -4.44033682e-05, +4.86234846e-06],
            [-1.02811799e-05, -4.06157099e-05, +6.02207637e-06],
        ]
    )

    d4 = disp.DFTD4Dispersion(mol, xc="r2SCAN")
    assert d4.kernel()[1] == approx(ref, abs=1.0e-7)


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_energy_hf():
    mol = gto.M(
        atom="""
             N  -1.57871857  -0.04661102   0.00000000
             N   1.57871857   0.04661102   0.00000000
             H  -2.15862174   0.13639605   0.80956529
             H  -0.84947130   0.65819321   0.00000000
             H  -2.15862174   0.13639605  -0.80956529
             H   2.15862174  -0.13639605  -0.80956529
             H   0.84947130  -0.65819321   0.00000000
             H   2.15862174  -0.13639605   0.80956529
             """
    )
    mf = disp.energy(scf.RHF(mol))
    assert mf.kernel() == approx(-110.91742452859162, abs=1.0e-8)


@pytest.mark.skipif(pyscf is None, reason="requires pyscf")
def test_gradient_hf():
    mol = gto.M(
        atom="""
             O  -1.65542061  -0.12330038   0.00000000
             O   1.24621244   0.10268870   0.00000000
             H  -0.70409026   0.03193167   0.00000000
             H  -2.03867273   0.75372294   0.00000000
             H   1.57598558  -0.38252146  -0.75856129
             H   1.57598558  -0.38252146   0.75856129
             """
    )
    ref = np.array(
        [
            [+1.72438133e-02, +5.08406920e-02, -1.89235207e-15],
            [+3.80018285e-02, -4.60223790e-02, +1.94297606e-16],
            [-3.05058266e-02, -1.26478132e-02, +1.04061735e-15],
            [+6.92338575e-03, -3.82898692e-02, +8.54978540e-16],
            [-1.58316004e-02, +2.30596847e-02, +2.18908543e-02],
            [-1.58316004e-02, +2.30596847e-02, -2.18908543e-02],
        ]
    )
    grad = disp.energy(scf.RHF(mol)).run().nuc_grad_method()
    assert grad.kernel() == approx(ref, abs=1.0e-8)
