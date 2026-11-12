import random

import pytest


@pytest.fixture
def reference():
    """Two deterministic random contigs. No data file is needed."""
    rng = random.Random(7)
    return {
        "chrA": "".join(rng.choice("ACGT") for _ in range(50_000)),
        "chrB": "".join(rng.choice("ACGT") for _ in range(30_000)),
    }
