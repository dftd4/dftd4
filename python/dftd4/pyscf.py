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
PySCF Support
-------------

Compatibility layer for supporting DFT-D4 in `pyscf <https://pyscf.org/>`_.
"""

try:
    from pyscf import lib, gto
except ModuleNotFoundError:
    raise ModuleNotFoundError("This submodule requires pyscf installed")

import numpy as np
from typing import Tuple

from .interface import DispersionModel, DampingParam


class DFTD4Dispersion(lib.StreamObject):
    """
    Implementation of the interface for using DFT-D4 in pyscf.

    Examples
    --------
    >>> from pyscf import gto
    >>> import dftd4.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          C   -0.755422531  -0.796459123  -1.023590391
    ...          C    0.634274834  -0.880017014  -1.075233285
    ...          C    1.406955202   0.199695367  -0.653144334
    ...          C    0.798863737   1.361204515  -0.180597909
    ...          C   -0.593166787   1.434312023  -0.133597923
    ...          C   -1.376239198   0.359205222  -0.553258516
    ...          I   -1.514344238   3.173268101   0.573601106
    ...          H    1.110906949  -1.778801728  -1.440619836
    ...          H    1.399172302   2.197767355   0.147412751
    ...          H    2.486417780   0.142466525  -0.689380574
    ...          H   -2.454252250   0.422581120  -0.512807958
    ...          H   -1.362353593  -1.630564523  -1.348743149
    ...          S   -3.112683203   6.289227834   1.226984439
    ...          H   -4.328789697   5.797771251   0.973373089
    ...          C   -2.689135032   6.703163830  -0.489062886
    ...          H   -1.684433029   7.115457372  -0.460265708
    ...          H   -2.683867206   5.816530502  -1.115183775
    ...          H   -3.365330613   7.451201412  -0.890098894
    ...          '''
    ... )
    >>> d4 = disp.DFTD4Dispersion(mol, xc="r2SCAN")
    >>> d4.kernel()[0]
    array(-0.0050011)
    """

    def __init__(self, mol, xc: str = "hf", atm: bool = True):
        self.mol = mol
        self.verbose = mol.verbose
        self.xc = xc
        self.atm = atm
        self.edisp = None
        self.grads = None

    def dump_flags(self, verbose=None) -> "DFTD4Dispersion":
        """
        Show options used for the DFT-D4 dispersion correction.
        """
        lib.logger.info(self, "** DFTD4 parameter **")
        lib.logger.info(self, "func %s", self.xc)
        return self

    def kernel(self) -> Tuple[float, np.ndarray]:
        """
        Compute the DFT-D4 dispersion correction.

        The dispersion model as well as the parameters are created locally and
        not part of the state of the instance.

        Returns
        -------
        float, ndarray
            The energy and gradient of the DFT-D4 dispersion correction.
        """
        mol = self.mol

        disp = DispersionModel(
            mol.atom_charges(),
            mol.atom_coords(),
            mol.charge,
        )

        param = DampingParam(
            method=self.xc,
            atm=self.atm,
        )

        res = disp.get_dispersion(param=param, grad=True)

        self.edisp = res.get("energy")
        self.grads = res.get("gradient")
        return self.edisp, self.grads

    def reset(self, mol) -> "DFTD4Dispersion":
        """
        Reset mol and clean up relevant attributes for scanner mode
        """
        self.mol = mol
        return self


class _DFTD4:
    """
    Stub class used to identify instances of the `DFTD4` class
    """

    pass


class _DFTD4Grad:
    """
    Stub class used to identify instances of the `DFTD4Grad` class
    """

    pass


def energy(mf):
    """
    Apply DFT-D4 corrections to SCF or MCSCF methods by returning an
    instance of a new class built from the original instances class.

    Parameters
    ----------
    mf
        The method to which DFT-D4 corrections will be applied.

    Returns
    -------
    The method with DFT-D4 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd4.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          N  -1.57871857  -0.04661102   0.00000000
    ...          N   1.57871857   0.04661102   0.00000000
    ...          H  -2.15862174   0.13639605   0.80956529
    ...          H  -0.84947130   0.65819321   0.00000000
    ...          H  -2.15862174   0.13639605  -0.80956529
    ...          H   2.15862174  -0.13639605  -0.80956529
    ...          H   0.84947130  -0.65819321   0.00000000
    ...          H   2.15862174  -0.13639605   0.80956529
    ...          '''
    ... )
    >>> mf = disp.energy(scf.RHF(mol)).run()
    converged SCF energy = -110.917424528592
    >>> mf.kernel()
    -110.917424528592
    """

    from pyscf.scf import hf
    from pyscf.mcscf import casci

    if not isinstance(mf, (hf.SCF, casci.CASCI)):
        raise TypeError("mf must be an instance of SCF or CASCI")

    with_dftd4 = DFTD4Dispersion(
        mf.mol,
        xc="hf"
        if isinstance(mf, casci.CASCI)
        else getattr(mf, "xc", "HF").upper().replace(" ", ""),
    )

    if isinstance(mf, _DFTD4):
        mf.with_dftd4 = with_dftd4
        return mf

    class DFTD4(_DFTD4, mf.__class__):
        """
        Patched SCF class including DFT-D4 corrections.
        """

        def __init__(self, method, with_dftd4: DFTD4Dispersion):
            self.__dict__.update(method.__dict__)
            self.with_dftd4 = with_dftd4
            self._keys.update(["with_dftd4"])

        def dump_flags(self, verbose=None) -> "DFTD4":
            mf.__class__.dump_flags(self, verbose)
            if self.with_dftd4:
                self.with_dftd4.dump_flags(verbose)
            return self

        def energy_nuc(self) -> float:
            enuc = mf.__class__.energy_nuc(self)
            if self.with_dftd4:
                enuc += self.with_dftd4.kernel()[0]
            return enuc

        def reset(self, mol=None) -> "DFTD4":
            self.with_dftd4.reset(mol)
            return mf.__class__.reset(self, mol)

        def nuc_grad_method(self) -> "DFTD4Grad":
            scf_grad = mf.__class__.nuc_grad_method(self)
            return grad(scf_grad)

        Gradients = lib.alias(nuc_grad_method, alias_name="Gradients")

    return DFTD4(mf, with_dftd4)


def grad(mfgrad):
    """
    Apply DFT-D4 corrections to SCF or MCSCF nuclear gradients methods
    by returning an instance of a new class built from the original class.

    Parameters
    ----------
    mfgrad
        The method to which DFT-D4 corrections will be applied.

    Returns
    -------
    The method with DFT-D4 corrections applied.

    Examples
    --------
    >>> from pyscf import gto, scf
    >>> import dftd4.pyscf as disp
    >>> mol = gto.M(
    ...     atom='''
    ...          O  -1.65542061  -0.12330038   0.00000000
    ...          O   1.24621244   0.10268870   0.00000000
    ...          H  -0.70409026   0.03193167   0.00000000
    ...          H  -2.03867273   0.75372294   0.00000000
    ...          H   1.57598558  -0.38252146  -0.75856129
    ...          H   1.57598558  -0.38252146   0.75856129
    ...          '''
    ... )
    >>> grad = disp.energy(scf.RHF(mol)).run().nuc_grad_method()
    converged SCF energy = -149.939098424774
    >>> g = grad.kernel()
    --------------- DFTD4 gradients ---------------
             x                y                z
    0 O     0.0172438133     0.0508406920     0.0000000000
    1 O     0.0380018285    -0.0460223790    -0.0000000000
    2 H    -0.0305058266    -0.0126478132    -0.0000000000
    3 H     0.0069233858    -0.0382898692    -0.0000000000
    4 H    -0.0158316004     0.0230596847     0.0218908543
    5 H    -0.0158316004     0.0230596847    -0.0218908543
    ----------------------------------------------
    """
    from pyscf.grad import rhf as rhf_grad

    if not isinstance(mfgrad, rhf_grad.Gradients):
        raise TypeError("mfgrad must be an instance of Gradients")

    # Ensure that the zeroth order results include DFTD4 corrections
    if not getattr(mfgrad.base, "with_dftd4", None):
        mfgrad.base = dftd4(mfgrad.base)

    class DFTD4Grad(_DFTD4Grad, mfgrad.__class__):
        """
        Patched SCF class including DFT-D4 corrections.
        """

        def grad_nuc(self, mol=None, atmlst=None):
            nuc_g = mfgrad.__class__.grad_nuc(self, mol, atmlst)
            with_dftd4 = getattr(self.base, "with_dftd4", None)
            if with_dftd4:
                disp_g = with_dftd4.kernel()[1]
                if atmlst is not None:
                    disp_g = disp_g[atmlst]
                nuc_g += disp_g
            return nuc_g

    dgrad = DFTD4Grad.__new__(DFTD4Grad)
    dgrad.__dict__.update(mfgrad.__dict__)
    return dgrad
