"""Focused Beam tests that exercise the custom-beam workflow.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import warnings

import numpy as np
import pytest
from scipy.constants import e as qe

try:
    from scipy.integrate import trapezoid
except ImportError:
    from scipy.integrate import trapz as trapezoid

    print(
        "[!] Using deprecated trapz function from scipy.integrate. Please update your scipy version."
    )

from bihc.beam import Beam


def _make_scheme(length, true_indices):
    scheme = [False] * length
    for idx in true_indices:
        scheme[idx] = True
    return scheme


@pytest.mark.filterwarnings("ignore:ppbk=.*:UserWarning")
def test_beam_custom_scheme_sets_profiles_and_charge():
    scheme = _make_scheme(6, [0, 2, 5])
    beam = Beam(
        M=len(scheme),
        fillingScheme=scheme,
        machine="PS",
        realMachineLength=False,
        spectrum="numeric",
        ppbk=12,
        bunchLength=1.0e-9,
        verbose=False,
    )

    assert beam.filledSlots == 3
    assert beam.totalBeamCharge == pytest.approx(beam.filledSlots * beam.Np * qe)

    t, profile = beam.longitudinalProfile
    assert t.shape == profile.shape
    area = trapezoid(beam.profile_1_bunch[1], beam.profile_1_bunch[0])
    assert area > 0

    samples_per_slot = len(beam.profile_1_bunch[0])
    for slot, filled in enumerate(scheme):
        start = slot * samples_per_slot
        segment = profile[start : start + samples_per_slot]
        if filled:
            np.testing.assert_array_equal(segment, beam.profile_1_bunch[1])
        else:
            assert not np.any(segment)


@pytest.mark.filterwarnings("ignore:ppbk=.*:UserWarning")
def test_beam_rejects_invalid_bunch_shape():
    scheme = _make_scheme(4, [0, 1, 2, 3])
    beam = Beam(
        M=len(scheme),
        fillingScheme=scheme,
        machine="PS",
        realMachineLength=False,
        spectrum="numeric",
        ppbk=8,
        bunchLength=1.0e-9,
        verbose=False,
    )

    with pytest.raises(ValueError):
        beam.bunchShape = "TRIANGLE"


def test_beam_rejects_removed_real_fft_option():
    with pytest.raises(ValueError, match="spectrum"):
        Beam(
            M=1,
            fillingScheme=[True],
            machine="PS",
            spectrum="numeric-rfft",
            verbose=False,
        )


def test_numeric_beam_warns_when_ppbk_undersamples_bunch_length():
    parameters = dict(
        M=2,
        fillingScheme=[True, True],
        machine="PS",
        realMachineLength=False,
        spectrum="numeric",
        bunchLength=np.array([40e-9, 20e-9]),
        verbose=False,
    )
    with pytest.warns(UserWarning, match="ppbk=16"):
        Beam(ppbk=16, **parameters)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        Beam(ppbk=24, **parameters)
    assert not caught


def test_numeric_beam_keeps_default_ppbk_without_explicit_fmax():
    beam = Beam(
        M=1,
        fillingScheme=[True],
        machine="PS",
        realMachineLength=False,
        t0=20e-9,
        bunchLength=4e-9,
        spectrum="numeric",
        verbose=False,
    )

    assert beam.ppbk == 250
    assert beam.spectrum[0][-1] > beam.fmax


@pytest.mark.parametrize(
    "ppbk, fmax, increased",
    [(None, 1e9, False), (8, 1e9, True), (1000, 1e9, False), (8, 0.1e9, True)],
)
def test_explicit_fmax_limits_numeric_spectrum_and_preserves_safe_ppbk(
    ppbk, fmax, increased
):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        beam = Beam(
            M=2,
            fillingScheme=[True, False],
            machine="PS",
            realMachineLength=False,
            t0=20e-9,
            bunchLength=4e-9,
            ppbk=ppbk,
            fmax=fmax,
            spectrum="numeric",
            verbose=False,
        )

    if ppbk is not None:
        assert beam.ppbk >= ppbk
    assert beam.ppbk >= 250
    assert any("ppbk increased" in str(w.message) for w in caught) == increased
    assert beam.t0 / beam.ppbk <= beam.sigma_time
    assert beam.ppbk / (2 * beam.t0) >= 1.5 * beam.fmax
    assert beam.spectrum[0][-1] <= beam.fmax * (1 + 1e-12)
    assert beam.fmax - beam.spectrum[0][-1] < beam.spectrum[0][1]


@pytest.mark.filterwarnings("ignore:ppbk=.*:UserWarning")
def test_set_bunches_overrides_profile_data():
    scheme = _make_scheme(4, [0])
    beam = Beam(
        M=len(scheme),
        fillingScheme=scheme,
        machine="PS",
        realMachineLength=False,
        spectrum="numeric",
        ppbk=4,
        bunchLength=1.0e-9,
        verbose=False,
    )

    t = np.linspace(0.0, 1e-6, beam.ppbk * 2)
    s = np.linspace(0.0, 1.0, beam.ppbk * 2)
    beam.setBunches([t, s], interp=False)

    assert np.allclose(beam.profile_1_bunch[0], t[: beam.ppbk])
    assert np.allclose(beam.profile_1_bunch[1], s[: beam.ppbk])
    assert np.allclose(beam.longitudinalProfile[0], t)
    assert np.allclose(beam.longitudinalProfile[1], s)
