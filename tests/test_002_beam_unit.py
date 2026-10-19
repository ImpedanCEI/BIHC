"""Focused Beam tests that exercise the custom-beam workflow.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import numpy as np
import pytest
from scipy.constants import e

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
    assert beam.totalBeamCharge == pytest.approx(beam.filledSlots * beam.Np * e)

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


def test_sampling_uses_shortest_bunch():
    parameters = dict(
        M=2,
        fillingScheme=[True, True],
        machine="PS",
        realMachineLength=False,
        spectrum="numeric",
        bunchLength=np.array([40e-9, 20e-9]),
        verbose=False,
    )
    Beam(ppbk=16, **parameters)
    Beam(ppbk=24, **parameters)


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
    assert beam.fmax == beam.ppbk / (2 * beam.t0)
    assert beam.spectrum[0][-1] <= beam.fmax


def test_fmax_adjusts_default_ppbk():
    beam = Beam(
        M=2,
        fillingScheme=[True, False],
        machine="PS",
        realMachineLength=False,
        t0=20e-9,
        bunchLength=4e-9,
        fmax=10e9,
        spectrum="numeric",
        verbose=False,
    )

    assert beam.ppbk == 400
    assert beam.ppbk / (2 * beam.t0) >= 10e9
    assert beam.fmax == beam.ppbk / (2 * beam.t0)
    assert beam.spectrum[0][-1] <= beam.fmax


def test_explicit_fmax_overrides_ppbk(capsys):
    beam = Beam(
        M=2,
        fillingScheme=[True, False],
        machine="PS",
        realMachineLength=False,
        t0=20e-9,
        bunchLength=4e-9,
        ppbk=25,
        fmax=1e9,
        spectrum="numeric",
        verbose=False,
    )
    output = capsys.readouterr().out

    assert beam.ppbk == 40
    assert "ppbk=25 overridden by fmax" in output
    assert beam.fmax == beam.ppbk / (2 * beam.t0)
    assert beam.spectrum[0][-1] > 0.5e9


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
