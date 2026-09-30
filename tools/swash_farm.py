#!/usr/bin/env python3
"""Run many SWASH instances in parallel for synthetic-dataset generation.

SWASH keeps all state in global module variables, so several instances can
only run as separate processes, each in its own working directory holding an
INPUT file.  This launcher

  1. renders an INPUT template (Python str.format placeholders such as {hs},
     {tp}, {seed}) once per parameter set read from a CSV/JSON file,
  2. generates the gridded input files with a user hook (default:
     tests/cases/gen_inputs.py <case> <dir>),
  3. runs N instances concurrently, pinning each to a GPU round-robin via
     CUDA_VISIBLE_DEVICES (start nvidia-cuda-mps-control -d first to share a
     GPU between several processes efficiently),
  4. writes runs/<name>/manifest.json with parameters, exit status and timing.

    swash_farm.py --exe build/bin/swash.exe --template tests/cases/marina_lay2/INPUT \\
                  --params params.csv --gen "python3 tests/cases/gen_inputs.py marina_lay2" \\
                  --jobs 4 --gpus 0,1 --out runs/dataset01

params.csv columns become template placeholders; a column named "id" names
the run directory (default: row index).
"""
import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor


def load_params(path):
    if path.endswith(".json"):
        return json.load(open(path))
    with open(path, newline="") as f:
        return [dict(r) for r in csv.DictReader(f)]


def run_one(i, params, args, gpu):
    name = str(params.get("id", f"run{i:04d}"))
    d = os.path.join(args.out, name)
    os.makedirs(d, exist_ok=True)
    with open(args.template) as f:
        text = f.read().format(**params)
    with open(os.path.join(d, "INPUT"), "w") as f:
        f.write(text)
    env = dict(os.environ)
    if gpu is not None:
        env["CUDA_VISIBLE_DEVICES"] = str(gpu)
    log = open(os.path.join(d, "stdout.log"), "w")
    t0 = time.time()
    rc = 0
    if args.gen:
        rc = subprocess.call(shlex.split(args.gen) + [d], stdout=log, stderr=subprocess.STDOUT, cwd=args.cwd)
    if rc == 0:
        rc = subprocess.call([os.path.abspath(args.exe)], cwd=d, env=env, stdout=log, stderr=subprocess.STDOUT)
    log.close()
    ok = rc == 0 and os.path.exists(os.path.join(d, "norm_end"))
    return {"id": name, "dir": d, "params": params, "gpu": gpu, "exit": rc, "ok": ok, "wall_s": round(time.time() - t0, 1)}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, help="path to swash.exe")
    ap.add_argument("--template", required=True, help="INPUT template with {placeholders}")
    ap.add_argument("--params", required=True, help="CSV or JSON list of parameter sets")
    ap.add_argument("--gen", default="", help="command run as '<gen> <rundir>' to create gridded inputs")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--gpus", default="", help="comma separated GPU ids for CUDA_VISIBLE_DEVICES round-robin")
    ap.add_argument("--out", default="runs")
    ap.add_argument("--cwd", default=os.getcwd())
    args = ap.parse_args()

    params = load_params(args.params)
    gpus = [g.strip() for g in args.gpus.split(",") if g.strip()] or [None]
    os.makedirs(args.out, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = [ex.submit(run_one, i, p, args, gpus[i % len(gpus)]) for i, p in enumerate(params)]
        for f in futs:
            r = f.result()
            results.append(r)
            print(f"{r['id']:<12s} gpu={r['gpu']} exit={r['exit']} ok={r['ok']} {r['wall_s']} s", flush=True)
    with open(os.path.join(args.out, "manifest.json"), "w") as f:
        json.dump({"exe": args.exe, "template": args.template, "runs": results}, f, indent=1)
    nfail = sum(not r["ok"] for r in results)
    print(f"== {len(results) - nfail}/{len(results)} runs ok; manifest in {args.out}/manifest.json")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
