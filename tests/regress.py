#!/usr/bin/env python3
"""Compare the outputs of a SWASH run against a golden (reference) run.

    regress.py <golden-dir> <run-dir> [--tier cpu|gpu] [--atol X] [--rtol Y]

Compared:
  * every field in fields.mat (MATLAB block output; variables such as
    Watlev_000000_000500, vel_x_..., vel_k1_...): max absolute error against
    atol + rtol * max|golden|, per quantity and over all times;
  * gauges.tbl (TABLE NOHEAD with TSEC WATL): the same tolerance on the
    time series, plus the significant wave height Hm0 and the peak period Tp
    from a Welch spectrum per gauge (physics-level check);
  * loads.tbl if present (floating-object forces).

Tiers (defaults, override with --atol/--rtol):
  cpu  : CPU-vs-CPU builds/compilers  atol 1e-4, rtol 1e-4, Hm0 0.5 %, Tp 1 %
  gpu  : GPU-vs-CPU "physics" tier    atol 1e-3, rtol 1e-3, Hm0 1 %,   Tp 2 %

Exit status 1 on any failure.  Requires numpy and scipy.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np

try:
    from scipy.io import loadmat
    from scipy.signal import welch
except ImportError:  # pragma: no cover
    sys.exit("regress.py needs scipy: pip install scipy")

TIERS = {
    "cpu": dict(atol=1e-4, rtol=1e-4, hm0=0.005, tp=0.01),
    "gpu": dict(atol=1e-3, rtol=1e-3, hm0=0.01, tp=0.02),
}
EXC = -99.0  # SWASH exception value for dry / undefined points


def load_fields(path):
    """Return {quantity: [(time_tag, array), ...]} from a SWASH .mat block file."""
    if not os.path.exists(path):
        return {}
    raw = loadmat(path)
    out = {}
    pat = re.compile(r"^(.+?)_(\d{6,8}_\d{3,6})$")
    for name, val in raw.items():
        if name.startswith("__"):
            continue
        m = pat.match(name)
        if m:
            q, tag = m.group(1), m.group(2)
        else:
            q, tag = name, "static"
        out.setdefault(q, []).append((tag, np.asarray(val, dtype=np.float64)))
    for q in out:
        out[q].sort(key=lambda t: t[0])
    return out


def load_table(path, ncols):
    """NOHEAD table: whitespace separated numbers, ncols per row."""
    if not os.path.exists(path):
        return None
    vals = np.loadtxt(path).reshape(-1, ncols)
    return vals


def compare_arrays(name, a, b, atol, rtol):
    a = np.where(a == EXC, np.nan, a)
    b = np.where(b == EXC, np.nan, b)
    if a.shape != b.shape:
        return False, f"{name}: shape mismatch {a.shape} vs {b.shape}"
    mask = ~(np.isnan(a) | np.isnan(b))
    if not mask.any():
        return True, f"{name}: no wet points to compare"
    diff = np.abs(a[mask] - b[mask]).max()
    scale = np.abs(b[mask]).max()
    tol = atol + rtol * scale
    ok = diff <= tol
    return ok, f"{name:<32s} max|diff|={diff:10.3e}  tol={tol:10.3e}  {'ok' if ok else 'FAIL'}"


def gauge_stats(t, eta):
    """Hm0 and Tp from a Welch spectrum of a surface elevation time series."""
    dt = np.median(np.diff(t))
    eta = eta - eta.mean()
    nper = min(len(eta), max(64, len(eta) // 4))
    f, S = welch(eta, fs=1.0 / dt, nperseg=nper)
    trap = getattr(np, "trapezoid", None) or np.trapz
    m0 = trap(S, f)
    hm0 = 4.0 * np.sqrt(max(m0, 0.0))
    fp = f[np.argmax(S[1:]) + 1] if len(f) > 1 else np.nan
    return hm0, (1.0 / fp if fp > 0 else np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("golden")
    ap.add_argument("run")
    ap.add_argument("--tier", choices=TIERS, default="cpu")
    ap.add_argument("--atol", type=float)
    ap.add_argument("--rtol", type=float)
    args = ap.parse_args()
    tol = dict(TIERS[args.tier])
    if args.atol is not None:
        tol["atol"] = args.atol
    if args.rtol is not None:
        tol["rtol"] = args.rtol

    failures = 0
    print(f"== regression {args.run} vs {args.golden} (tier {args.tier}: {tol})")

    # ---- gridded fields
    g = load_fields(os.path.join(args.golden, "fields.mat"))
    r = load_fields(os.path.join(args.run, "fields.mat"))
    for q in sorted(g):
        if q not in r:
            print(f"{q}: missing in run"); failures += 1; continue
        gl, rl = g[q], r[q]
        if len(gl) != len(rl):
            print(f"{q}: {len(gl)} frames in golden vs {len(rl)} in run"); failures += 1
        worst = (0.0, "")
        for (tg, ag), (tr, ar) in zip(gl, rl):
            ok, msg = compare_arrays(f"{q}@{tg}", ag, ar, tol["atol"], tol["rtol"])
            if not ok:
                failures += 1; print(msg)
            d = float(re.search(r"max\|diff\|=\s*([0-9.e+-]+)", msg).group(1)) if "max|diff|" in msg else 0.0
            if d >= worst[0]:
                worst = (d, msg)
        print(f"  {q:<28s} frames={len(gl):4d}  worst: {worst[1].split(':',1)[-1].strip() if worst[1] else 'n/a'}")

    # ---- gauges
    gt = load_table(os.path.join(args.golden, "gauges.tbl"), 2)
    rt = load_table(os.path.join(args.run, "gauges.tbl"), 2)
    if gt is not None and rt is not None:
        # rows are (time, location) blocks; infer number of gauges from repeated times
        times = gt[:, 0]
        ngauge = int(np.sum(times == times[0]))
        gt = gt.reshape(-1, ngauge, 2); rt = rt.reshape(-1, ngauge, 2)
        ok, msg = compare_arrays("gauges WATL", gt[:, :, 1], rt[:, :, 1], tol["atol"], tol["rtol"])
        print("  " + msg)
        if not ok:
            failures += 1
        for i in range(ngauge):
            hg, tg_ = gauge_stats(gt[:, i, 0], gt[:, i, 1])
            hr, tr_ = gauge_stats(rt[:, i, 0], rt[:, i, 1])
            eh = abs(hr - hg) / max(hg, 1e-6)
            et = abs(tr_ - tg_) / max(tg_, 1e-6) if np.isfinite(tg_) and np.isfinite(tr_) else 0.0
            flag = "ok" if (eh <= tol["hm0"] and et <= tol["tp"]) else "FAIL"
            if flag == "FAIL":
                failures += 1
            print(f"  gauge {i+1}: Hm0 {hg:.4f} -> {hr:.4f} ({100*eh:.2f} %)   Tp {tg_:.3f} -> {tr_:.3f} ({100*et:.2f} %)  {flag}")
    elif gt is not None:
        print("  gauges.tbl missing in run"); failures += 1

    # ---- floating-object loads (optional)
    gl = load_table(os.path.join(args.golden, "loads.tbl"), 7)
    rl = load_table(os.path.join(args.run, "loads.tbl"), 7)
    if gl is not None and rl is not None:
        ok, msg = compare_arrays("loads", gl[:, 1:], rl[:, 1:], tol["atol"] * 1e3, tol["rtol"])
        print("  " + msg)
        if not ok:
            failures += 1

    print(f"== {'PASS' if failures == 0 else f'FAIL ({failures} checks)'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
