#!/usr/bin/env python3
"""Rebuild the immune-context figure layouts from deposited source data."""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

FIGURES = Path(__file__).resolve().parent
sys.path.insert(0, str(FIGURES.parents[1] / "release"))
from layout import replay_view


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="Extracted v1.2.0 processed-data deposit")
    parser.add_argument("output", type=Path, help="New or empty historical-figure output directory")
    args = parser.parse_args()
    data, output = args.data.resolve(), args.output.resolve()
    if output.is_relative_to(data) or output.is_relative_to(FIGURES.parents[2]):
        parser.error("Keep output outside the source and data directories")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        parser.error("Choose a new or empty output directory")
    with replay_view(data) as compatibility:
        render(compatibility, output)


def render(data, output):
    results = data / "data/r10/results"
    if not results.is_dir():
        raise ValueError("Missing the immune-context results in the replay view")
    output.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix="asd-context-figure-replay-") as work:
        stage = Path(work)
        shutil.copytree(results, stage / "03_results")
        environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "MPLCONFIGDIR": str(stage / "mpl")}
        for script in ("render_context_main.py", "render_context_supplement.py", "context_figure_documentation.py"):
            subprocess.run([sys.executable, FIGURES / script, "--data", data, "--stage", stage],
                           env=environment, check=True)
        for folder in ("04_figures", "05_supplementary"):
            shutil.move(stage / folder, output / folder)


if __name__ == "__main__":
    main()
