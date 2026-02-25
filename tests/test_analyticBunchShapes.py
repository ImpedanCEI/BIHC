"""
Test for analytic spectrum with
different profile shape

* date: 23/03/2026
* author: Elena de la Fuente
"""

from pathlib import Path

import pytest

import bihc

_EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
CSV_FILE = (
    _EXAMPLES
    / "BunchProfilesComparison"
    / "25ns_2760b_2748_2494_2572_288bpi_13inj.csv"
)
IMPEDANCE_FILE = _EXAMPLES / "MinMaxDissipatedPower" / "PillboxImpedance.txt"


def test_analytic_gaussian():
    prof = "GAUSSIAN"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        verbose=False,
        spectrum="analytic",
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(2129.9978925335763, rel=1e-3)


def test_analytic_parabolic():
    prof = "PARABOLIC"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        verbose=False,
        spectrum="analytic",
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(4764.592327462002, rel=1e-3)


def test_analytic_cos2():
    prof = "COS2"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        verbose=False,
        spectrum="analytic",
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(3276.6520910812233, rel=1e-3)
