# -*- coding: utf-8 -*-
"""Common I/O utilities for single-cell analysis pipeline."""
import os
from pathlib import Path
from datetime import datetime


def ensure_dir(path):
    """Create directory if it doesn't exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def get_output_path(base_dir, filename):
    """Get full output path, creating directory if needed."""
    ensure_dir(base_dir)
    return os.path.join(base_dir, filename)


def log_io(operation, input_path, output_path, shape_before=None, shape_after=None):
    """Log I/O operation details."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    msg = f"[{timestamp}] {operation}\n"
    msg += f"  Input: {input_path}\n"
    msg += f"  Output: {output_path}\n"
    if shape_before:
        msg += f"  Shape before: {shape_before}\n"
    if shape_after:
        msg += f"  Shape after: {shape_after}\n"
    return msg


def check_file_exists(path):
    """Check if file exists, raise error if not."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"Required file not found: {path}")
    return True


def get_file_size_kb(path):
    """Get file size in KB."""
    return os.path.getsize(path) / 1024