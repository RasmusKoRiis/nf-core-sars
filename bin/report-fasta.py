#!/usr/bin/env python3
"""Compatibility entry point for FASTA-only reports.

The FASTA workflow now uses the same merge/report logic as the FASTQ workflow
so resistance flags and QC-derived fields stay consistent between modes.
"""
from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    runpy.run_path(str(Path(__file__).with_name("report.py")), run_name="__main__")
