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

"""
Tests for the ASE interface to dftd4.
"""

def test_d4_molecular():
    """Short test for DFT-D4 on a subset of the G2."""
    from pytest import approx
    from ase.collections import g2
    from dftd4.calculators import D4_model
    thr = 1.0e-7
    subset = {'C4H4NH': (-0.17254765346459802, 0.0016481180998958495),
              'CH3SCH3': (-0.118396636155385, 0.0013823525532589362),
              'SiH2_s3B1d': (-0.018605813985880824, 0.00015462713987174245),
              'CH3SH': (-0.05722345063877332, 0.0005448462666493717),
              'CH3CO': (-0.05123732582276016, 0.0007742597318805521),
              'C3H4_D2d': (-0.07030306476081483, 0.0009590568552901135),
              'ClF3': (-0.04213170121261114, 0.0007279615294069727),
              'H2CO': (-0.01994226769935472, 0.0002148654202479407),
              'CH3COOH': (-0.08745304310824802, 0.00114445197261909),
              'HCF3': (-0.03887439301151592, 0.00029142286765719783),
              'CH3S': (-0.045072896650802385, 0.0003779814622251629),
              'C3H6_D3h': (-0.09235218055900413, 0.0006339294824490826),
              'CS2': (-0.06617104839166844, 0.001825513384064356),
              'SiH2_s1A1d': (-0.018788257002315355, 0.0003824353060351632),
              'C4H4S': (-0.19935202085994883, 0.0018416124706281756)}

    calc = D4_model(s6=1.0, s8=1.26829475, a1=0.39907098, a2=5.03951304, s9=1.0)
    for name, data in subset.items():
        atoms = g2[name]
        energy, gnorm = data
        atoms.set_calculator(calc)
        assert atoms.get_potential_energy() == approx(energy, thr)
        assert abs(atoms.get_forces().flatten()).mean() == approx(gnorm, thr)

def test_d3_molecular():
    """Short test for DFT-D4 on a subset of the G2."""
    from pytest import approx
    from ase.collections import g2
    from dftd4.calculators import D3_model
    thr = 1.0e-7
    subset = {'cyclobutene': (-0.1458649601010507, 0.0038664473279371447),
              'CH3ONO': (-0.07478103530920273, 0.0013917101416313475),
              'SiH3': (-0.03286579996845799, 0.0005983134328171301),
              'C3H6_D3h': (-0.09950955015994396, 0.0013365866217120727),
              'CO2': (-0.021752475167461482, 0.00019725833839493446),
              'ClF3': (-0.035917624197226, 0.0005044951380930323),
              'C3H4_D2d': (-0.078068075574534, 0.0013471319373085805),
              'COF2': (-0.029919370204709785, 0.00027306557116953074),
              '2-butyne': (-0.12279742624908908, 0.0019018023476061548),
              'C2H5': (-0.0583205245708783, 0.0008123602763700175),
              'BF3': (-0.02553364750247437, 0.00026111693285554887),
              'SiH2_s3B1d': (-0.021066607316984308, 0.0003363117490443568),
              'N2O': (-0.022068066394571734, 0.0002284681148899124)}

    calc = D3_model(s6=1.0, s8=1.44956635, a1=0.32704579, a2=5.36988013, s9=1.0)
    for name, data in subset.items():
        atoms = g2[name]
        energy, gnorm = data
        atoms.set_calculator(calc)
        assert atoms.get_potential_energy() == approx(energy, thr)
        assert abs(atoms.get_forces().flatten()).mean() == approx(gnorm, thr)
