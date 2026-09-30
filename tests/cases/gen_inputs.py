#!/usr/bin/env python3
"""Generate the gridded input files (bathymetry, draft) for the reference cases.

The files are small enough to regenerate on every run, so they are not committed.
Layout: IDLA=3 as used by the READINP commands in the decks, i.e. the first line
is the row y = 0 and x varies fastest along a line.  Bottom values are depths
below the still water level, positive downward (SWASH convention).

Usage: gen_inputs.py <case_name> <output_dir>
"""
import sys
import numpy as np


def write_field(path, field):
    # field[j, i] with j the y index, i the x index
    with open(path, "w") as f:
        for row in field:
            f.write(" ".join(f"{v:.4f}" for v in row) + "\n")


def marina_bathymetry(mx=200, my=150, dx=1.0):
    """6 m deep at the west boundary, shoaling to 2.5 m in the east half,
    a dry quay in the north-east corner and a breakwater strip.  Values on
    the (mx+1) x (my+1) input grid points."""
    x = np.arange(mx + 1) * dx
    y = np.arange(my + 1) * dx
    X, Y = np.meshgrid(x, y)
    depth = np.where(X < 80.0, 6.0, 6.0 - 3.5 * np.clip((X - 80.0) / 60.0, 0.0, 1.0))
    # quay: dry land (bottom above still water level -> negative depth)
    quay = (X > 150.0) & (Y > 110.0)
    depth[quay] = -2.0
    # detached breakwater: emergent strip
    bw = (X > 120.0) & (X < 124.0) & (Y > 40.0) & (Y < 90.0)
    depth[bw] = -1.5
    return depth


def basin_bathymetry(mx=300, my=300, dx=2.0):
    """Flat 8 m basin with a gentle 1:50 beach in the east third."""
    x = np.arange(mx + 1) * dx
    y = np.arange(my + 1) * dx
    X, Y = np.meshgrid(x, y)
    depth = 8.0 - np.clip((X - 400.0) / 50.0, 0.0, None)
    return depth


def pontoon_draft(mx=200, my=150, dx=1.0, x0=60.0, y0=60.0, L=30.0, B=10.0, draft=1.5):
    """Draft field for a fixed rectangular pontoon.  Cells outside the object
    get a bottom-of-object level 10 m ABOVE the still water level (value -10),
    which never pressurises the flow (see SwashPresFlow: pressurised if
    s1 >= -flos)."""
    x = np.arange(mx + 1) * dx
    y = np.arange(my + 1) * dx
    X, Y = np.meshgrid(x, y)
    inside = (X >= x0) & (X <= x0 + L) & (Y >= y0) & (Y <= y0 + B)
    return np.where(inside, draft, -10.0)


def main():
    case, out = sys.argv[1], sys.argv[2]
    if case in ("marina_lay2", "marina_float"):
        write_field(f"{out}/bot.dat", marina_bathymetry())
        if case == "marina_float":
            write_field(f"{out}/draft.dat", pontoon_draft())
    elif case == "basin_da":
        write_field(f"{out}/bot.dat", basin_bathymetry())
    else:
        sys.exit(f"unknown case {case}")


if __name__ == "__main__":
    main()
