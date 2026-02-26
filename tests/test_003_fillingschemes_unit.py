"""Unit tests for the public filling scheme helpers.

The scenarios mirror the hand-crafted schemes used in the examples
folder, but keep the inputs intentionally small so the behaviour can
be asserted deterministically.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import pytest

from bihc.fillingschemes import (
    fillingSchemeLHC,
    fillingSchemeLHC_8b4e,
    fillingSchemeLHC_standard,
    fillingSchemeSPS,
    fillingSchemeSPS_8b4e,
    fillingSchemeSPS_AWAKE,
    fillingSchemeSPS_BCMS,
    fillingSchemeSPS_standard,
)


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


def test_sps_standard_total_slots():
    ntrains = 3
    scheme = fillingSchemeSPS_standard(ntrains)

    assert len(scheme) == 920
    assert sum(scheme) == ntrains * 72
    assert not any(scheme[-20:])


def test_sps_bcms_keeps_batch_gaps():
    scheme = fillingSchemeSPS_BCMS(ntrains=2)

    assert len(scheme) == 920
    assert scheme[:48] == [True] * 48
    assert scheme[48:56] == [False] * 8


def test_sps_8b4e_pattern_repeats_correctly():
    scheme = fillingSchemeSPS_8b4e(ntrains=1)
    expected_cell = [True] * 8 + [False] * 4
    assert scheme[:12] == expected_cell


def test_sps_awake_single_bunch():
    scheme = fillingSchemeSPS_AWAKE(nbunches=1, ntrains=1)
    assert scheme[0] is True
    assert sum(scheme) == 1


def test_lhc_standard_counts_true_slots():
    ninj, nbunches, ntrains = 2, 5, 1
    scheme = fillingSchemeLHC_standard(
        ninj=ninj,
        nbunches=nbunches,
        ntrains=ntrains,
        batchspacing=3,
        injspacing=5,
    )

    assert len(scheme) == 3564
    assert sum(scheme) == ninj * ntrains * nbunches


def test_lhc_8b4e_injects_expected_pattern():
    scheme = fillingSchemeLHC_8b4e(ninj=1, ntrains=1)
    cell = [True] * 8 + [False] * 4
    assert scheme[:12] == cell
    assert len(scheme) == 3564
