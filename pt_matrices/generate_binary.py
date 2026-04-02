#!/usr/bin/env python
"""Generate binary (.bin) versions of all PT matrix text files.

For each .dat file loaded by CLASS-PT, this script reads the whitespace-
separated doubles and writes them as a raw double-precision binary file
with suffix .dat.bin (matching the LOAD_MATRIX macro convention).

Usage:
    python generate_binary.py [--force]

Options:
    --force   Overwrite existing .bin files (default: skip if exists)
"""
import glob
import os
import sys

import numpy as np

CLASSDIR = os.path.dirname(os.path.abspath(__file__))
ROOTDIR = os.path.dirname(CLASSDIR)

force = "--force" in sys.argv

# Directories containing matrix .dat files
search_dirs = [
    os.path.join(ROOTDIR, "pt_matrices"),
    os.path.join(ROOTDIR, "pt_matrices", "compute_matrices_python"),
    os.path.join(ROOTDIR, "ORTHOGONAL_NG_matrices", "computation_of_matrices"),
]

count = 0
skipped = 0

for search_dir in search_dirs:
    if not os.path.isdir(search_dir):
        continue
    for datfile in sorted(glob.glob(os.path.join(search_dir, "**", "*.dat"), recursive=True)):
        # Skip files that are themselves .bin files
        if datfile.endswith(".bin"):
            continue
        binfile = datfile + ".bin"
        if os.path.exists(binfile) and not force:
            skipped += 1
            continue
        try:
            data = np.loadtxt(datfile)
            data = data.ravel().astype(np.float64)
            data.tofile(binfile)
            count += 1
            print(f"  wrote {binfile} ({len(data)} doubles, {os.path.getsize(binfile)} bytes)")
        except Exception as e:
            print(f"  SKIP {datfile}: {e}")

print(f"\nDone: {count} files written, {skipped} skipped (already exist)")
