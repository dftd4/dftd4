#!/usr/bin/env python3

import subprocess, sys, json, os, pytest

if len(sys.argv) < 5:
    raise RuntimeError("Requires at least four arguments")

prog = sys.argv[1]
geom = sys.argv[2]
outp = sys.argv[3]
args = sys.argv[4:]

stat = subprocess.call(
    [prog, geom, "--json", os.path.basename(outp)] + args,
    shell=False,
    stdin=None,
    stderr=subprocess.STDOUT,
)
if stat != 0:
    raise RuntimeError("Calculation failed")

with open(outp) as f:
    ref = json.load(f)
    del ref["version"]

with open(os.path.basename(outp)) as f:
    res = json.load(f)

for key in ref:
    if key not in res:
        raise RuntimeError("Missing '" + key + "' entry in results")
    assert pytest.approx(res[key], abs=1.0e-9) == ref[key]
