import os
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
BAMS = [str(DATA / "nov20_barcode09_subset.bam"), str(DATA / "nov20_barcode10_subset.bam")]

LAB_MODKIT = "/usr/xtmp/nd141/programs/dist_modkit_v0.4.4_7cf558c/modkit"


def find_modkit():
    for c in (os.environ.get("FIBERNHAT_MODKIT"), shutil.which("modkit"), LAB_MODKIT):
        if c and os.path.exists(c):
            return c
    return None


@pytest.fixture(scope="session")
def bams():
    return BAMS


@pytest.fixture(scope="session")
def modkit():
    mk = find_modkit()
    if mk is None:
        import warnings
        warnings.warn("modkit not found: PARITY TESTS SKIPPED. Install ont-modkit=0.4.4 or set FIBERNHAT_MODKIT.")
        pytest.skip("modkit not installed")
    return mk
