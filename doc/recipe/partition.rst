How to distribute a calculation over MPI ranks?
===============================================

When every MPI rank holds the same structure, each rank would normally repeat
the complete dispersion calculation. A work partition assigns a disjoint share
of the pairwise and ATM interaction loops to each rank instead.

Parts are zero based and every unit of work belongs to exactly one part, so
summing the energy, gradient and virial over all parts reproduces the complete
result. DFT-D4 itself performs no communication; the reduction is left to the
caller.

.. note::

   Structure-dependent quantities such as coordination numbers, charges and
   :math:`C_6` coefficients are evaluated for the full system on every part.
   The speedup is therefore bounded by the interaction loops, which dominate
   for larger systems and when the ATM contribution is enabled.

Partition the interaction loops
-------------------------------

The Fortran API takes the partition as an optional argument; omitting it selects
the complete work. The C API stores it on the dispersion model, next to the
cutoffs.

.. tab-set::
   :sync-group: code

   .. tab-item:: Fortran
      :sync: fortran

      .. literalinclude:: partition-example/mpi.f90
         :language: fortran
         :caption: mpi.f90

   .. tab-item:: C
      :sync: c

      .. literalinclude:: partition-example/mpi.c
         :language: c
         :caption: mpi.c

The examples can be compiled using ``mpifort`` or ``mpicc`` together with the
flags reported by ``pkg-config dftd4 mctc-lib --cflags --libs``. The result is
independent of the number of ranks up to the summation order.

Running with one part is identical to omitting the partition. In Fortran this
is also available as ``serial_work_partition``; in C and Python the model can
be reset with ``part=0`` and ``nparts=1``.

The numerical Hessian uses the same partition. The pairwise decomposition and
the model properties do not.
