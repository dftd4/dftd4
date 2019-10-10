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

"""Tests for the ctypes interface to dftd4."""


def test_library():
    """check if we can find the library and it looks okay"""
    from ctypes import cdll
    from ctypes.util import find_library

    name = find_library("dftd4")
    assert name is not None
    # check if we can load this one
    lib = cdll.LoadLibrary(name)
    # check if we find some functions
    assert lib.D4_calculation is not None
    assert lib.D4_PBC_calculation is not None
    assert lib.D3_calculation is not None
    assert lib.D3_PBC_calculation is not None

    from dftd4.interface import DFTD4Library

    libdftd4 = DFTD4Library(library=lib)
    assert libdftd4.library is lib
    # check if we find some functions
    assert hasattr(libdftd4, "D4Calculation")
    assert hasattr(libdftd4, "D3Calculation")

    libdftd4 = DFTD4Library()
    # check if we find some functions
    assert hasattr(libdftd4, "D4Calculation")
    assert hasattr(libdftd4, "D3Calculation")


class Testing:
    """test all defined interfaces"""
    from ctypes import c_int, c_double
    from numpy import array
    from dftd4.interface import DFTD4Library

    natoms = 24
    numbers = array(
        [6, 7, 6, 7, 6, 6, 6, 8, 7, 6, 8, 7, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        dtype=c_int,
    )
    positions = array(
        [
            [2.02799738646442, 0.09231312124713, -0.14310895950963],
            [4.75011007621000, 0.02373496014051, -0.14324124033844],
            [6.33434307654413, 2.07098865582721, -0.14235306905930],
            [8.72860718071825, 1.38002919517619, -0.14265542523943],
            [8.65318821103610, -1.19324866489847, -0.14231527453678],
            [6.23857175648671, -2.08353643730276, -0.14218299370797],
            [5.63266886875962, -4.69950321056008, -0.13940509630299],
            [3.44931709749015, -5.48092386085491, -0.14318454855466],
            [7.77508917214346, -6.24427872938674, -0.13107140408805],
            [10.30229550927022, -5.39739796609292, -0.13672168520430],
            [12.07410272485492, -6.91573621641911, -0.13666499342053],
            [10.70038521493902, -2.79078533715849, -0.14148379504141],
            [13.24597858727017, -1.76969072232377, -0.14218299370797],
            [7.40891694074004, -8.95905928176407, -0.11636933482904],
            [1.38702118184179, 2.05575746325296, -0.14178615122154],
            [1.34622199478497, -0.86356704498496, 1.55590600570783],
            [1.34624089204623, -0.86133716815647, -1.84340893849267],
            [5.65596919189118, 4.00172183859480, -0.14131371969009],
            [14.67430918222276, -3.26230980007732, -0.14344911021228],
            [13.50897177220290, -0.60815166181684, 1.54898960808727],
            [13.50780014200488, -0.60614855212345, -1.83214617078268],
            [5.41408424778406, -9.49239668625902, -0.11022772492007],
            [8.31919801555568, -9.74947502841788, 1.56539243085954],
            [8.31511620712388, -9.76854236502758, -1.79108242206824],
        ],
        dtype=c_double,
    )
    lib = DFTD4Library()

    def test_d4_interface(self):
        """check if the D4 interface is working correctly."""
        thr = 1.0e-8
        from pytest import approx
        from dftd4.interface import P_MBD_APPROX_ATM, P_REFQ_EEQ

        options = {
            'lmbd': P_MBD_APPROX_ATM,
            'refq': P_REFQ_EEQ,
            'wf': 6.0,
            'g_a': 3.0,
            'g_c': 2.0,
            'properties': True,
            'energy': True,
            'forces': True,
            'hessian': False,
            'print_level': 1,
            's6': 1.0000,  # B3LYP-D4-ATM parameters
            's8': 1.93077774,
            's9': 1.0,
            's10': 0.0,
            'a1': 0.40520781,
            'a2': 4.46255249,
            'alp': 16,
        }
        kwargs = {
            "natoms": self.natoms,
            "numbers": self.numbers,
            "charge": 0.0,
            "positions": self.positions,
            "options": options,
            "output": "-",
        }
        results = self.lib.D4Calculation(**kwargs)
        assert results
        gnorm = abs(results["gradient"].flatten()).mean()
        assert approx(results["energy"], thr) == -0.04928342326381067
        assert approx(gnorm, thr) == 0.0002453003088176903

    def test_d3_interface(self):
        """check if the D3 interface is working correctly."""
        thr = 1.0e-8
        from pytest import approx
        from dftd4.interface import P_MBD_APPROX_ATM, P_REFQ_EEQ

        options = {
            'lmbd': P_MBD_APPROX_ATM,
            'refq': P_REFQ_EEQ,
            'wf': 4.0,
            'g_a': 0.0,
            'g_c': 0.0,
            'properties': True,
            'energy': True,
            'forces': True,
            'hessian': False,
            'print_level': 1,
            's6': 1.0000,  # B3LYP-D3-ATM parameters
            's8': 2.17707910,
            's9': 1.0,
            's10': 0.0,
            'a1': 0.31888091,
            'a2': 4.85097930,
            'alp': 16,
        }
        kwargs = {
            "natoms": self.natoms,
            "numbers": self.numbers,
            "positions": self.positions,
            "options": options,
            "output": "-",
        }
        results = self.lib.D3Calculation(**kwargs)
        assert results
        gnorm = abs(results["gradient"].flatten()).mean()
        assert approx(results["energy"], thr) == -0.049193282449077516
        assert approx(gnorm, thr) == 0.0003855527267206884


def test_load_D4_damping_parameters():
    """check if we can load parameters from the shared object"""
    from dftd4.interface import DFTD4Library

    funcs = {'pbe0': {'s6': 1.0, 's8': 1.20065498, 's10': 0.0, 'a1': 0.40085597,
                      'a2': 5.02928789, 's9': 1.0, 'alp': 16, 'beta': 1.0},
             'b3lyp': {'s6': 1.0, 's8': 2.02929367, 's10': 0.0, 'a1': 0.40868035,
                       'a2': 4.53807137, 's9': 1.0, 'alp': 16, 'beta': 1.0},
             'B3-LYP': {'s6': 1.0, 's8': 2.02929367, 's10': 0.0, 'a1': 0.40868035,
                        'a2': 4.53807137, 's9': 1.0, 'alp': 16, 'beta': 1.0},
             'PwpB95': {'s6': 0.82, 's8': -0.34639127, 's10': 0.0,
                        'a1': 0.41080636, 'a2': 3.83878274, 's9': 1.0, 'alp': 16,
                        'beta': 1.0},
             'Lh07sSVWN': {'s6': 1.0, 's8': 3.16675531, 's10': 0.0,
                           'a1': 0.35965552, 'a2': 4.31947614, 's9': 1.0,
                           'alp': 16, 'beta': 1.0},
             'revPBE': {'s6': 1.0, 's8': 1.7467653, 's10': 0.0, 'a1': 0.536349,
                        'a2': 3.07261485, 's9': 1.0, 'alp': 16, 'beta': 1.0}}

    library = DFTD4Library()
    for func in funcs:
        assert funcs[func] == library.D4Parameters(func)


def test_structs():
    """test if structures are robust against input"""
    from dftd4.interface import DFTDParameters, DFTDOptions

    parameters = {'s6': 1.0, 's8': 1.44956635, 'a1': 0.32704579, 'a2': 5.36988013,
                  's9': 1.0, 'alp': 16, 's10': 0.0, 'beta': 1.0,
                  'wrong_key': 'stuff'}
    dparam = DFTDParameters(**parameters).to_dict()
    del parameters['wrong_key']
    assert parameters == dparam

    parameters = {'s6': 1.0, 's8': 1.44956635, 'a1': 0.32704579, 'a2': 5.36988013,
                  's9': 1.0, 'alp': 16}
    dparam = DFTDParameters(**parameters).to_dict()
    assert parameters == {key: dparam[key] for key in parameters}

    options = {'lmbd': 3, 'refq': 5, 'wf': 6.0, 'g_a': 3.0, 'g_c': 2.0,
               'properties': True, 'energy': True, 'forces': True,
               'hessian': False, 'print_level': 2,
               'wrong_key': 'stuff'}

    dopt = DFTDOptions(**options).to_dict()
    del options['wrong_key']
    assert options == dopt

    options = {'lmbd': 3, 'refq': 5, 'wf': 6.0, 'g_a': 3.0, 'g_c': 2.0,
               'properties': True, 'energy': True, 'forces': True}

    dopt = DFTDOptions(**options).to_dict()
    assert options == {key: dopt[key] for key in options}
