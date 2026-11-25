"""Resolve inputs from the extracted processed-data deposit."""
from pathlib import Path


def released_file(data_dir, relative_path):
    root = Path(data_dir).resolve(strict=True)
    candidate = root / "data" / "r7" / relative_path
    try:
        resolved = candidate.resolve(strict=True)
    except FileNotFoundError as error:
        raise FileNotFoundError(
            f"Required release input missing: data/r7/{relative_path}"
        ) from error
    if not resolved.is_relative_to(root) or not resolved.is_file():
        raise ValueError(f"Input must be a release file: data/r7/{relative_path}")
    return resolved
