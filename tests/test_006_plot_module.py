"""Minimal tests for the plotting utilities.

These tests do not assert on pixels; they merely ensure the public
APIs execute without raising and produce artefacts such as files or
progress text when expected.

* date: 25/02/2026
* author: Elena de la Fuente
"""

import io

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np
import pytest

from bihc.impedance import Impedance
from bihc.plot import Plot, progressbar


class PlotHarness(Plot):
    def __init__(self):
        self.longitudinalProfile = [
            np.linspace(0.0, 1e-6, 5),
            np.linspace(1.0, 5.0, 5),
        ]
        self.filledSlots = 3
        self._fillNumber = 7
        self.fmax = 4e9
        self._spectrum = [
            np.linspace(0.0, 4e9, 5),
            np.linspace(0.1, 0.5, 5),
        ]

    @property
    def spectrum(self):
        return self._spectrum


@pytest.fixture
def plot_harness():
    return PlotHarness()


@pytest.fixture
def show_spy(monkeypatch):
    calls = []

    def fake_show(*args, **kwargs):
        calls.append((args, kwargs))

    monkeypatch.setattr(plt, "show", fake_show)
    return calls


def test_progressbar_prints_completion():
    stream = io.StringIO()
    list(progressbar(range(3), prefix="test", size=3, out=stream))
    assert "100%" in stream.getvalue()


def test_plot_longitudinal_profile_invokes_show(plot_harness, show_spy):
    plot_harness.plotLongitudinalProfile()
    assert len(show_spy) == 1


def test_plot_power_spectrum_adds_label(plot_harness, show_spy):
    plot_harness._fillNumber = 123
    plot_harness.plotPowerSpectrum(fmin=0, fmax=2)
    assert len(show_spy) == 1


def test_plot_impedance_draws_curves(show_spy):
    freqs = np.linspace(0.0, 1e9, 5)
    Z = Impedance(f=freqs)
    Z.Zr = np.linspace(0.0, 1.0, 5)
    Z.Zi = -Z.Zr
    Z.plotImpedance()
    assert len(show_spy) == 1


def test_plot_spectrum_and_impedance_handles_interpolation(
    plot_harness, show_spy
):
    freqs = np.linspace(0.0, 5e9, 6)
    Z = Impedance(f=freqs)
    Z.Zr = np.linspace(0.2, 1.2, 6)
    Z.Zi = np.zeros_like(Z.Zr)
    plot_harness.plotSpectrumAndImpedance(Z)
    assert len(show_spy) == 1


def test_save_longitudinal_distribution_writes_file(plot_harness, tmp_path):
    target = tmp_path / "profile.txt"
    plot_harness.saveLongitudinalDistribution(
        target, phase_shift=5, normalise=True
    )
    data = np.loadtxt(target, skiprows=1)
    assert data.shape[1] == 2
    # Normalisation ensures the maximum is 1
    assert np.max(data[:, 1]) == pytest.approx(1.0)


def test_plot2beam_handles_shift(plot_harness, show_spy):
    partner = PlotHarness()
    Plot.plot2Beam(
        plot_harness, partner, shift=plot_harness.longitudinalProfile[0][1]
    )
    assert len(show_spy) == 1
