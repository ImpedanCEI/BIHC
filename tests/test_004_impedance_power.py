"""Unit coverage for impedance helpers and the power mixin.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import numpy as np
import pytest
from scipy.constants import c
from scipy.constants import e as qe

from bihc import power as power_module
from bihc.beam import Beam
from bihc.impedance import Impedance


def _beam_with_custom_spectrum(
    freqs, amps, time_span, filled_slots=1, np_bunch=1e11
):
    beam = Beam(
        M=1,
        fillingScheme=[True],
        machine="PS",
        realMachineLength=False,
        spectrum="user",
        ppbk=4,
        bunchLength=1e-9,
        verbose=False,
    )
    beam.longitudinalProfile = [
        np.array([0.0, time_span], dtype=float),
        np.zeros(2, dtype=float),
    ]
    beam.setSpectrum(
        [np.array(freqs, dtype=float), np.array(amps, dtype=float)]
    )
    beam.filledSlots = filled_slots
    beam.Np = np_bunch
    return beam


def _impedance_from(freqs, real=None):
    Z = Impedance(f=freqs)
    if real is None:
        real = np.ones_like(freqs)
    Z.Zr = np.array(real, dtype=float)
    Z.Zi = np.zeros_like(Z.Zr)
    return Z


@pytest.fixture
def noop_tqdm(monkeypatch):
    monkeypatch.setattr(
        power_module, "tqdm", lambda iterable, *_, **__: iterable
    )


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
    beam = _beam_with_custom_spectrum(
        freqs=freqs, amps=amps, time_span=1e-6, filled_slots=2, np_bunch=1e11
    )

    Z = Impedance(f=freqs)
    Z.Zr = np.array([0.0, 50.0])
    Z.Zi = np.zeros_like(Z.Zr)

    loss, density = beam.getPloss(Z)

    mask = freqs > 0.0
    S = amps[mask]
    Zr = Z.Zr[mask]
    f0 = 1 / (beam.longitudinalProfile[0][-1] - beam.longitudinalProfile[0][0])
    expected_density = 2 * (f0 * qe * S * beam.filledSlots * beam.Np) ** 2 * Zr

    assert loss == pytest.approx(expected_density.sum())
    assert density[~mask] == pytest.approx(0.0)
    assert np.allclose(density[mask], expected_density)


def test_shifted_ploss_scans_frequency_grid(noop_tqdm):
    freqs = np.linspace(0.0, 3e6, 4)
    amps = np.linspace(0.1, 0.4, 4)
    beam = _beam_with_custom_spectrum(
        freqs=freqs, amps=amps, time_span=1e-6, filled_slots=1, np_bunch=1e10
    )

    Z = _impedance_from(freqs, real=np.linspace(1.0, 2.0, len(freqs)))
    shifts, power = beam.getShiftedPloss(Z, shift=freqs[1] - freqs[0])

    # size = 1 -> shifts covers [-1, 0]
    assert np.array_equal(shifts, np.array([-1, 0]))
    assert power.shape == shifts.shape
    assert hasattr(beam, "Zmax")
    assert np.array_equal(beam.Zmax.f, freqs)


def test_shifted_power_spectrum_returns_density_grid(noop_tqdm):
    freqs = np.linspace(0.0, 4e6, 5)
    amps = np.linspace(0.2, 1.0, 5)
    beam = _beam_with_custom_spectrum(
        freqs=freqs, amps=amps, time_span=2e-6, filled_slots=2, np_bunch=5e9
    )

    Z = _impedance_from(freqs, real=np.linspace(0.5, 2.5, len(freqs)))
    shifts, power_spectrum = beam.getShiftedPowerSpectrum(
        Z, shift=freqs[1] - freqs[0]
    )

    assert power_spectrum.shape == (len(shifts), len(freqs))
    assert np.all(power_spectrum >= 0)


def test_two_beam_ploss_accepts_tau_list(noop_tqdm):
    freqs = np.array([0.0, 2.0])
    amps = np.array([0.0, 1.0])
    beam = _beam_with_custom_spectrum(
        freqs=freqs, amps=amps, time_span=1.0, filled_slots=1, np_bunch=1e9
    )

    Z0 = _impedance_from(freqs, real=np.array([0.0, 10.0]))
    tau = np.array([0.0, 1e-9])

    results = beam.get2BeamPloss(Z0, tau_s=tau)

    assert results[0] == pytest.approx(0.0)
    assert results[1] > results[0]
    assert hasattr(beam, "P_loss")


def test_two_beam_ploss_with_offsets_and_secondary_impedance(noop_tqdm):
    freqs = np.array([0.0, 2.0])
    amps = np.array([0.0, 0.5])
    beam = _beam_with_custom_spectrum(
        freqs=freqs, amps=amps, time_span=1.0, filled_slots=2, np_bunch=5e9
    )

    Z0 = _impedance_from(freqs, real=np.array([0.0, 5.0]))
    Z1 = _impedance_from(freqs, real=np.array([0.0, 2.0]))
    s = np.array([0.0, 10.0])  # meters
    offset1 = np.array([0.0, 0.002])
    offset2 = np.array([0.001, 0.0])

    tau = 2 * s / c
    results = beam.get2BeamPloss(
        Z0,
        s=s,
        tau_s=tau,
        offset1=offset1,
        offset2=offset2,
        Z_1=Z1,
    )

    assert len(results) == len(s)
    assert results[-1] > results[0]
