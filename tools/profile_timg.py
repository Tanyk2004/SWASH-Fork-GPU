#!/usr/bin/env python3
"""Parse the timing tables that a TIMG build writes at the end of the PRINT file.

Two tables are printed:
  * the aggregate table from SWPRTI ("# Details on timings of the simulation")
  * the per-section table from SwashPrintTimers ("# Per-section timers"), which
    lists every timer id that accumulated time, with the section name, the
    number of calls and the cpu / wall-clock seconds.

    profile_timg.py PRINT [--top N]
"""
import re
import sys

agg_re = re.compile(r"^\s*\d+\s+#\s+(.*?):\s+([0-9.]+)\s+([0-9.]+)\s*$")
sec_re = re.compile(r"^\s*\d+\s+#\s+(\d+)\s+(\S.*?)\s{2,}(\d+)\s+([0-9.]+)\s+([0-9.]+)\s*$")


def main():
    path = sys.argv[1]
    top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 25
    agg, sec = [], []
    for line in open(path, errors="replace"):
        m = agg_re.match(line)
        if m:
            agg.append((m.group(1).strip(), float(m.group(2)), float(m.group(3))))
            continue
        m = sec_re.match(line)
        if m:
            sec.append((int(m.group(1)), m.group(2).strip(), int(m.group(3)), float(m.group(4)), float(m.group(5))))
    if agg:
        total = next((w for n, c, w in agg if n.startswith("total time")), None)
        print(f"{'aggregate section':<30s} {'cpu [s]':>10s} {'wall [s]':>10s} {'% wall':>7s}")
        for n, c, w in agg:
            pct = f"{100*w/total:6.1f}" if total else ""
            print(f"{n:<30s} {c:10.2f} {w:10.2f} {pct:>7s}")
    if sec:
        total = max(w for *_, w in sec)
        print()
        print(f"{'id':>4s} {'section':<36s} {'calls':>8s} {'cpu [s]':>10s} {'wall [s]':>10s} {'% wall':>7s}")
        for i, n, k, c, w in sorted(sec, key=lambda t: -t[4])[:top]:
            print(f"{i:4d} {n:<36s} {k:8d} {c:10.2f} {w:10.2f} {100*w/total:7.1f}")
    if not agg and not sec:
        print("no timing table found (build with -DTIMG=ON)")


if __name__ == "__main__":
    main()
