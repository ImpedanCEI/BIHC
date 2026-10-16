"""Unit tests for the public filling scheme helpers.

The scenarios mirror the hand-crafted schemes used in the examples
folder, but keep the inputs intentionally small so the behaviour can
be asserted deterministically.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import pytest

from bihc.fillingschemes import (
    fillingSchemeFCCee,
    fillingSchemeFCCee_ttbar,
    fillingSchemeFCCee_WW,
    fillingSchemeFCCee_ZH,
    fillingSchemeLHC,
    fillingSchemeLHC_8b4e,
    fillingSchemeLHC_standard,
    fillingSchemeSPS,
    fillingSchemeSPS_8b4e,
    fillingSchemeSPS_AWAKE,
    fillingSchemeSPS_BCMS,
    fillingSchemeSPS_SFTPRO,
    fillingSchemeSPS_standard,
)


def test_fccee_scheme_keeps_trains_and_gaps():
    scheme = fillingSchemeFCCee(
        nslots=22,
        ntrains=2,
        nbunches=3,
        batchspacing=2,
        ninj=2,
        injspacing=1,
    )

    injection = [True] * 3 + [False] * 2 + [True] * 3 + [False]
    assert scheme == injection * 2 + [False] * 4


def test_fccee_scheme_rejects_oversized_pattern():
    with pytest.raises(ValueError, match="available machine slots"):
        fillingSchemeFCCee(nslots=10, ntrains=2, nbunches=5, batchspacing=1)


def test_fccee_scheme_defaults_to_study_pattern():
    scheme = fillingSchemeFCCee()

    assert len(scheme) == 12120
    assert sum(scheme) == 40 * 280
    assert scheme[:280] == [True] * 280
    assert scheme[280:303] == [False] * 23


@pytest.mark.parametrize(
    (
        "scheme_factory",
        "ntrains",
        "groups_per_train",
        "group_size",
        "group_gap",
        "train_gap",
    ),
    [
        (fillingSchemeFCCee_ZH, 2, 75, 2, 19, 3764),
        (fillingSchemeFCCee_WW, 8, 58, 4, 16, 371),
        (fillingSchemeFCCee_ttbar, 2, 16, 2, 151, 3763),
    ],
)
def test_fccee_mode_patterns(
    scheme_factory, ntrains, groups_per_train, group_size, group_gap, train_gap
):
    scheme = scheme_factory()
    train_length = groups_per_train * group_size + (groups_per_train - 1) * group_gap
    next_train = train_length + train_gap

    assert len(scheme) == 12120
    assert sum(scheme) == ntrains * groups_per_train * group_size
    assert scheme[:group_size] == [True] * group_size
    assert scheme[group_size : group_size + group_gap] == [False] * group_gap
    assert (
        scheme[group_size + group_gap : 2 * group_size + group_gap]
        == [True] * group_size
    )
    assert scheme[train_length:next_train] == [False] * train_gap
    assert scheme[next_train : next_train + group_size] == [True] * group_size


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
    second_train = scheme[nbunches + batchspacing : 2 * nbunches + batchspacing]

    assert first_train == [True] * nbunches
    assert gap == [False] * batchspacing
    assert second_train == [True] * nbunches


def test_sps_scheme_defaults_to_lhc25ns_study_pattern():
    scheme = fillingSchemeSPS()

    assert len(scheme) == 924
    assert sum(scheme) == 4 * 72
    assert scheme[72:81] == [False] * 9


def test_sps_standard_total_slots():
    ntrains = 3
    scheme = fillingSchemeSPS_standard(ntrains)

    assert len(scheme) == 924
    assert sum(scheme) == ntrains * 72
    assert not any(scheme[-20:])
    assert scheme == fillingSchemeSPS(ntrains=ntrains)


def test_sps_standard_respects_bunches_per_train():
    scheme = fillingSchemeSPS_standard(ntrains=2, nbunches=36)

    assert len(scheme) == 924
    assert sum(scheme) == 72
    assert scheme == fillingSchemeSPS(nbunches=36, ntrains=2)


def test_sps_bcms_keeps_batch_gaps():
    scheme = fillingSchemeSPS_BCMS(ntrains=2)

    assert len(scheme) == 924
    assert scheme[:48] == [True] * 48
    assert scheme[48:56] == [False] * 8
    assert scheme == fillingSchemeSPS(nbunches=48, ntrains=2, batchspacing=8)


def test_sps_bcms_respects_bunches_per_train():
    scheme = fillingSchemeSPS_BCMS(ntrains=2, nbunches=36)

    assert len(scheme) == 924
    assert sum(scheme) == 72
    assert scheme == fillingSchemeSPS(nbunches=36, ntrains=2, batchspacing=8)


def test_sps_8b4e_pattern_repeats_correctly():
    scheme = fillingSchemeSPS_8b4e(ntrains=1)
    expected_cell = [True] * 8 + [False] * 4
    assert scheme[:12] == expected_cell
    assert len(scheme) == 924
    assert sum(scheme) == 56
    assert not any(scheme[-20:])


def test_sps_awake_single_bunch():
    scheme = fillingSchemeSPS_AWAKE(nbunches=1, ntrains=1)
    assert scheme[0] is True
    assert sum(scheme) == 1


def test_sps_sftpro_thresholds_measured_intensities():
    intensities = [0.0] * 4620
    intensities[0] = 1e9
    intensities[1] = 1.1e9
    scheme = fillingSchemeSPS_SFTPRO(intensities)

    assert len(scheme) == 4620
    assert sum(scheme) == 1
    assert scheme[1] is True
    with pytest.raises(ValueError, match="4620 slots"):
        fillingSchemeSPS_SFTPRO(intensities[:-1])


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
    assert scheme == fillingSchemeLHC(
        ninj=ninj,
        nbunches=nbunches,
        ntrains=ntrains,
        batchspacing=3,
        injspacing=5,
    )


def test_lhc_8b4e_injects_expected_pattern():
    scheme = fillingSchemeLHC_8b4e(ninj=1, ntrains=1)
    cell = [True] * 8 + [False] * 4
    assert scheme[:12] == cell
    assert len(scheme) == 3564
