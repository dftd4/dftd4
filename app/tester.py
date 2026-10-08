#!/usr/bin/env python3
"""
Minimal Python wrapper for testing the dftd4 command line interface.

The wrapper will assume a specific order in the arguments rather than
providing a generic command line interface by itself since it is
supposed to be used by meson for testing purposes only.
"""

try:
    import subprocess
    import sys
    import json
    import os
    import pytest
except ImportError:
    raise SystemExit(77)

if len(sys.argv) < 4:
    raise RuntimeError("Requires at least four arguments")

thr = 1.0e-9
prog = sys.argv[1]
outp = sys.argv[2]
args = sys.argv[3:]

run = subprocess.run(
    [prog, "--json", os.path.basename(outp)] + args,
    shell=False,
    stdin=None,
    stderr=subprocess.STDOUT,
    stdout=subprocess.PIPE,
    universal_newlines=True,
)
print(run.stdout, end="")
if run.returncode != 0:
    raise RuntimeError("Calculation failed")

verbosity = (
    2
    + sum(arg in ("-v", "--verbose") for arg in args)
    - sum(arg in ("-s", "--silent") for arg in args)
)
assert ("Timing (wall time):" in run.stdout) == (verbosity > 0)
if "--hessian" in args and verbosity > 1:
    assert "numerical two-body Hessian" in run.stdout
    assert "analytical ATM Hessian" in run.stdout
    assert "numerical full Hessian" not in run.stdout

with open(outp) as f:
    ref = json.load(f)
    del ref["version"]

with open(os.path.basename(outp)) as f:
    res = json.load(f)

for key in ref:
    if key not in res:
        raise RuntimeError("Missing '" + key + "' entry in results")
    assert pytest.approx(res[key], abs=thr) == ref[key]
