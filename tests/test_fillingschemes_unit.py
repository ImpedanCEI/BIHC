"""Unit tests for the public filling scheme helpers.

The scenarios mirror the hand-crafted schemes used in the examples
folder, but keep the inputs intentionally small so the behaviour can
be asserted deterministically.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import pytest

from bihc.fillingschemes import fillingSchemeLHC, fillingSchemeSPS


def test_lhc_scheme_counts_true_slots():
    """The LHC helper should fill ninj * ntrains * nbunches buckets."""
    ninj, nbunches, ntrains = 2, 3, 2
    scheme = fillingSchemeLHC(
        ninj=ninj,
        nbunches=nbunches,
        ntrains=ntrains,
        batchspacing=1,
        injspacing=2,
    )

    assert len(scheme) == 3564
    assert sum(scheme) == ninj * ntrains * nbunches
    # The tail should be left as padding False slots.
    assert not any(scheme[-50:])


def test_lhc_scheme_raises_when_pattern_exceeds_ring():
    """When the pattern is longer than the machine, the helper raises."""
    with pytest.raises(Exception) as excinfo:
        fillingSchemeLHC(
            ninj=60,
            nbunches=200,
            ntrains=4,
            batchspacing=1,
            injspacing=1,
        )

    assert "available machine slots" in str(excinfo.value)


def test_sps_scheme_preserves_train_structure():
    """SPS helper keeps each train contiguous with batch gaps in between."""
    nbunches, ntrains, batchspacing = 5, 2, 3
    scheme = fillingSchemeSPS(
        nbunches=nbunches, ntrains=ntrains, batchspacing=batchspacing
    )

    first_train = scheme[:nbunches]
    gap = scheme[nbunches : nbunches + batchspacing]
    second_train = scheme[
        nbunches + batchspacing : 2 * nbunches + batchspacing
    ]

    assert first_train == [True] * nbunches
    assert gap == [False] * batchspacing
    assert second_train == [True] * nbunches
