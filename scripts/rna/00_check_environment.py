#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Step 0: Check environment
Verifies Python path, package versions, working directories, and data files.
"""
import sys
from pathlib import Path

def main():
    print("=" * 50)
    print("ENVIRONMENT CHECK")
    print("=" * 50)

    # Python path
    print(f"\nPython executable: {sys.executable}")

    # Package versions
    packages = ['scanpy', 'anndata', 'numpy', 'pandas', 'scipy', 'matplotlib']
    for pkg in packages:
        try:
            mod = __import__(pkg)
            ver = getattr(mod, '__version__', 'unknown')
            print(f"  {pkg}: {ver}")
        except ImportError:
            print(f"  {pkg}: NOT INSTALLED")

    # Working directory
    cwd = Path.cwd()
    print(f"\nWorking directory: {cwd}")

    # Check key directories
    dirs_to_check = [
        'data/raw/rna', 'data/processed/rna',
        'scripts/rna', 'scripts/common',
        'results/figures/rna', 'results/tables/rna',
        'configs', 'docs'
    ]
    print("\nDirectory structure:")
    for d in dirs_to_check:
        p = Path(d)
        status = "OK" if p.exists() else "MISSING"
        print(f"  {d}: {status}")

    # Check input data
    data_file = Path("data/raw/rna/pbmc3k_raw.h5ad")
    print(f"\nInput data: {data_file}")
    if data_file.exists():
        size_kb = data_file.stat().st_size / 1024
        print(f"  Size: {size_kb:.1f} KB - OK")
    else:
        print("  MISSING!")

    print("\n" + "=" * 50)
    print("CHECK COMPLETE")
    print("=" * 50)

if __name__ == "__main__":
    main()