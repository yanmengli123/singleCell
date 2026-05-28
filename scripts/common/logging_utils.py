# -*- coding: utf-8 -*-
"""Logging utilities for single-cell analysis pipeline."""
import os
import sys
from datetime import datetime
from pathlib import Path


def setup_logging(log_dir, log_name):
    """Set up logging to file and console."""
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    log_path = os.path.join(log_dir, log_name)

    # Create logger
    logger = Logger(log_path)

    # Also print to console
    logger.log_to_console = True

    return logger


class Logger:
    """Simple logger for pipeline steps."""

    def __init__(self, log_path):
        self.log_path = log_path
        self.log_to_console = False
        self.start_time = datetime.now()

        # Write header
        with open(self.log_path, 'w') as f:
            f.write(f"=== Pipeline Log ===\n")
            f.write(f"Start: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Python: {sys.executable}\n")
            f.write(f"Working dir: {os.getcwd()}\n")
            f.write("=" * 50 + "\n\n")

    def log(self, message):
        """Log a message."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        line = f"[{timestamp}] {message}"

        with open(self.log_path, 'a') as f:
            f.write(line + "\n")

        if self.log_to_console:
            print(line)

    def log_step(self, step_name, **kwargs):
        """Log a pipeline step with key-value pairs."""
        self.log(f"\n--- {step_name} ---")
        for key, value in kwargs.items():
            self.log(f"  {key}: {value}")

    def finish(self):
        """Log completion."""
        end_time = datetime.now()
        duration = end_time - self.start_time

        with open(self.log_path, 'a') as f:
            f.write("\n" + "=" * 50 + "\n")
            f.write(f"Finish: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Duration: {duration}\n")

        self.log(f"\nCompleted in {duration}")