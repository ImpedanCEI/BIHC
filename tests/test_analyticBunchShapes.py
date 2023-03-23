'''
Test for analytic spectrum with 
different profile shape

* date: 23/03/2023
* author: Elena de la Fuente
'''

import sys
sys.path.append('../../')

import bihc


def test_analytic_gaussian():
    prof = 'GAUSSIAN'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'

    beam = bihc.Beam(LPCfile=file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Get spectrum and profile
    [f, S] = beam_analytic.spectrum
    [t, s] = beam_analytic.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(impedance_file)

    # Compute power loss
    power = beam_analytic.getPloss(Z)[0]

    assert power == 2181.14217410196

def test_analytic_parabolic():
    prof = 'PARABOLIC'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'

    beam = bihc.Beam(LPCfile=file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Get spectrum and profile
    [f, S] = beam_analytic.spectrum
    [t, s] = beam_analytic.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(impedance_file)

    # Compute power loss
    power = beam_analytic.getPloss(Z)[0]

    assert power == 4878.996917588338

def test_analytic_parabolic():
    prof = 'PARABOLIC'
    file='25ns_2760b_2748_2494_2572_288bpi_13inj.csv'
    impedance_file = 'PillboxImpedance.txt'

    beam = bihc.Beam(LPCfile=file, bunchShape=prof, verbose=False, spectrum='analytic')

    # Get spectrum and profile
    [f, S] = beam_analytic.spectrum
    [t, s] = beam_analytic.profile_1_bunch

    # Importing an impedance curve
    Z = bihc.Impedance(f) 
    Z.getImpedanceFromCST(impedance_file)

    # Compute power loss
    power = beam_analytic.getPloss(Z)[0]

    assert power == 3355.3292665992517


