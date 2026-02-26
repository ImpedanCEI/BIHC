"""
Test for numeric spectrum with
different profile shape

* date: 23/03/2023
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


def test_numeric_gaussian():

    prof = "GAUSSIAN"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        ppbk=1000,
        verbose=False,
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(2131.886538626302, rel=1e-3)


def test_numeric_qgaussian():
    prof = "q-GAUSSIAN"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        ppbk=1000,
        verbose=False,
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(2773.0557710660837, rel=1e-3)


def test_numeric_parabolic():
    prof = "PARABOLIC"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        ppbk=1000,
        verbose=False,
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(1372.5726614514037, rel=1e-3)


def test_numeric_cos2():
    prof = "COS2"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        ppbk=1000,
        verbose=False,
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(1614.1388070946766, rel=1e-3)


def test_numeric_binomial():
    prof = "BINOMIAL"
    beam = bihc.Beam(
        LPCfile=str(CSV_FILE),
        bunchShape=prof,
        ppbk=1000,
        verbose=False,
    )

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f)
    Z.getImpedanceFromCST(str(IMPEDANCE_FILE))

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(1618.5694534921295, rel=1e-3)
