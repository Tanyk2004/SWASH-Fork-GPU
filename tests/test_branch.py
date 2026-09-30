#!/usr/bin/env python3
"""Checks that the in-memory checkpoint of the engine is complete:

    run A:  init, step N1, save slot 1, step N2      -> fields A
    run B:  restore slot 1, step N2                  -> fields B

A and B must be bitwise identical (the restored state must contain everything
the time stepping depends on).  Any difference means an array or scalar is
missing from engine/swash_state.ftn90 (regenerate with engine/gen_state.py
after adding it to the module lists) or that some routine keeps hidden state.

    test_branch.py <libswash_engine.so> <case> [--n1 200] [--n2 100]

Runs in tests/work/branch-<case>/ using the reference case's INPUT.
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(here, "..", "engine", "python"))
from swash_engine import Swash  # noqa: E402

FIELDS = ["WATL", "U", "V", "W", "Q", "HK", "WETS", "PRESP"]


def snapshot(sw):
    out = {}
    for f in FIELDS:
        for k in range(0, sw.nk + 1):
            try:
                out[f"{f}{k}"] = sw.field(f, k)
            except Exception:
                pass
    out["t"] = np.array([sw.time, sw.dt])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("lib")
    ap.add_argument("case")
    ap.add_argument("--n1", type=int, default=200)
    ap.add_argument("--n2", type=int, default=100)
    a = ap.parse_args()

    work = os.path.join(here, "work", f"branch-{a.case}")
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    shutil.copy(os.path.join(here, "cases", a.case, "INPUT"), os.path.join(work, "INPUT"))
    subprocess.check_call([sys.executable, os.path.join(here, "cases", "gen_inputs.py"), a.case, work])

    sw = Swash(os.path.abspath(a.lib))
    sw.init(work)
    assert sw.step(a.n1) == a.n1, "run ended before n1 steps"
    sw.save(1)
    assert sw.step(a.n2) == a.n2
    A = snapshot(sw)
    sw.restore(1)
    assert sw.step(a.n2) == a.n2
    B = snapshot(sw)
    sw.finish()

    bad = 0
    for k in A:
        if not np.array_equal(A[k], B[k]):
            d = np.nanmax(np.abs(A[k].astype(np.float64) - B[k].astype(np.float64)))
            print(f"MISMATCH {k}: max|diff| = {d:.3e}")
            bad += 1
    print("== branch test", "PASS" if bad == 0 else f"FAIL ({bad} fields differ)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
