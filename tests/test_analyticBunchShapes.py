'''
Test for analytic spectrum with 
different profile shape

* date: 23/03/2023
* author: Elena de la Fuente
'''

import bihc
import pytest

def test_analytic_gaussian():
    prof = 'GAUSSIAN'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == pytest.approx(2181.14217410196, 0.1) #TODO use approx for all the asserts

def test_analytic_parabolic():
    prof = 'PARABOLIC'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 4878.996917588338

def test_analytic_cos2():
    prof = 'COS2'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'
    path = 'tests/'
    beam = bihc.Beam(LPCfile=path+file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Get spectrum and profile
    [f, S] = beam.spectrum
    [t, s] = beam.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(path+impedance_file)

    # Compute power loss
    power = beam.getPloss(Z)[0]

    assert power == 3355.3292665992517


