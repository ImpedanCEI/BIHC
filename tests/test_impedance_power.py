"""Unit coverage for impedance helpers and the power mixin.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import numpy as np
import pytest
from scipy.constants import e as qe

from bihc.impedance import Impedance
from bihc.power import Power


class PowerHarness(Power):
    """Minimal Power subclass that exposes a programmable spectrum."""

    def __init__(self, freqs, amps, time_span, filled_slots=1, np_bunch=1e11):
        self.longitudinalProfile = [np.array([0.0, time_span]), np.zeros(2)]
        self._spectrum = [np.array(freqs), np.array(amps)]
        self.filledSlots = filled_slots
        self.Np = np_bunch
        self.verbose = False

    @property
    def spectrum(self):
        return self._spectrum


def test_resonator_impedance_real_part_at_resonance():
    freqs = np.array([1e8, 2e8, 3e8])
    Z = Impedance(f=freqs)
    Z.getResonatorImpedance(Rs=5.0, Qr=50.0, fr=2e8)

    center = np.argmin(np.abs(freqs - 2e8))
    assert Z.Zr[center] == pytest.approx(5.0)
    assert Z.Zi[center] == pytest.approx(0.0, abs=1e-12)


def test_impedance_addition_interpolates_to_reference_grid():
    ref = Impedance(f=np.linspace(1e6, 4e6, 4))
    ref.getResonatorImpedance(Rs=10.0, Qr=5.0, fr=2e6)

    coarse = Impedance(f=np.linspace(1e6, 4e6, 2))
    coarse.getRWImpedance(L=1.0, b=0.02, sigma=5.7e7)

    combined = ref + coarse
    expected = ref.Zr + np.interp(ref.f, coarse.f, coarse.Zr)

    assert np.allclose(combined.Zr, expected)


def test_power_loss_matches_manual_formula():
    freqs = np.array([0.0, 5e6])
    amps = np.array([0.0, 2.0])
    harness = PowerHarness(
        freqs=freqs, amps=amps, time_span=1e-6, filled_slots=2, np_bunch=1e11
    )

    Z = Impedance(f=freqs)
    Z.Zr = np.array([0.0, 50.0])
    Z.Zi = np.zeros_like(Z.Zr)

    loss, density = harness.getPloss(Z)

    mask = freqs > 0.0
    S = amps[mask]
    Zr = Z.Zr[mask]
    f0 = 1 / (
        harness.longitudinalProfile[0][-1] - harness.longitudinalProfile[0][0]
    )
    expected_density = (
        2 * (f0 * qe * S * harness.filledSlots * harness.Np) ** 2 * Zr
    )

    assert loss == pytest.approx(expected_density.sum())
    assert density[~mask] == pytest.approx(0.0)
    assert np.allclose(density[mask], expected_density)
