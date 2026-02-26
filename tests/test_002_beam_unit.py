"""Focused Beam tests that exercise the custom-beam workflow.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import numpy as np
import pytest
from scipy.constants import e as qe

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
    assert beam.totalBeamCharge == pytest.approx(
        beam.filledSlots * beam.Np * qe
    )

    t, profile = beam.longitudinalProfile
    assert t.shape == profile.shape
    area = np.trapezoid(beam.profile_1_bunch[1], beam.profile_1_bunch[0])
    assert area > 0


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
