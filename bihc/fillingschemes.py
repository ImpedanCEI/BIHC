# copyright ################################# #
# This file is part of the BIHC Package.      #
# Copyright (c) CERN, 2024.                   #
# ########################################### #

"""Filling schemes as boolean slot sequences.

True indicates a bunch; False indicates an empty slot.
"""

from __future__ import annotations

from collections.abc import Sequence

# FCC-ee filling scheme


def fillingSchemeFCCee(
    nslots: int = 12120,
    ntrains: int = 40,
    nbunches: int = 280,
    batchspacing: int = 23,
    ninj: int = 1,
    injspacing: int = 0,
) -> list[bool]:
    """Return the FCC-ee BPM-study train pattern by default.

    The default has 40 trains of 280 bunches on a 12,120-slot ring, with
    23 empty slots between trains. Only bunch occupancy is represented.

    Parameters
    ----------
    nslots : int, default 12120
        Total number of slots in the ring.
    ntrains : int, default 40
        Number of trains per injection.
    nbunches : int, default 280
        Number of filled slots per train.
    batchspacing : int, default 23
        Number of empty slots between trains.
    ninj : int, default 1
        Number of injections.
    injspacing : int, default 0
        Number of empty slots after each injection.

    Returns
    -------
    list of bool
        Occupied and empty slots across the ring.
    """
    if nslots <= 0 or ntrains <= 0 or nbunches <= 0 or ninj <= 0:
        raise ValueError("nslots, ntrains, nbunches, and ninj must be positive")
    if batchspacing < 0 or injspacing < 0:
        raise ValueError("batchspacing and injspacing must be nonnegative")

    train = [True] * nbunches
    gap = [False] * batchspacing
    injection = (train + gap) * (ntrains - 1) + train + [False] * injspacing
    filling = injection * ninj
    if len(filling) > nslots:
        raise ValueError(
            f"Filling scheme length: {len(filling)} > available machine slots: {nslots}"
        )
    return filling + [False] * (nslots - len(filling))


# LHC filling schemes


def fillingSchemeLHC(
    ninj: int, nbunches: int, ntrains: int, batchspacing: int = 7, injspacing: int = 37
) -> list[bool]:
    """Return an LHC 25 ns filling with configurable trains and injections.

    Each injection contains ``ntrains`` trains of ``nbunches`` filled slots.
    This is the general LHC bunch-train pattern used by the standard helper.

    Parameters
    ----------
    ninj : int
        Number of injections.
    nbunches : int
        Number of filled slots per train.
    ntrains : int
        Number of trains per injection.
    batchspacing : int, default 7
        Number of empty slots between trains.
    injspacing : int, default 37
        Number of empty slots after each injection.

    Returns
    -------
    list of bool
        Filling of the 3,564-slot LHC ring.
    """

    # Define filling scheme: parameters
    # ninj = 14 # Defining number of injections
    nslots = 3564  # Defining total number of slots for LHC
    # ntrain = 4 # Defining the number of trains
    # nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72..
    batchS = batchspacing  # Batch spacing in 25 ns slots
    injspacing = injspacing  # Injection spacing in 25 ns slots

    # Defining the trains as lists of True/Falses
    bt = [True] * nbunches
    st = [False] * batchS
    stt = [False] * injspacing
    sc = [False] * (
        nslots
        - (
            ntrains * nbunches * ninj
            + ((ntrains - 1) * (batchS) * ninj)
            + ((1) * injspacing * (ninj))
        )
    )
    an1 = (bt + st) * (ntrains - 1) + bt + stt
    an = (
        an1 * ninj + sc
    )  # This is the final true false sequence that is the beam distribution

    if len(an) > nslots:
        raise Exception(
            f"Filling scheme length: {len(an)} > available machine slots: {nslots}"
        )

    return an


def fillingSchemeLHC_standard(
    ninj: int,
    nbunches: int = 72,
    ntrains: int = 4,
    batchspacing: int = 7,
    injspacing: int = 37,
) -> list[bool]:
    """Return the standard LHC 25 ns multi-injection train pattern.

    Each injection contains four 72-bunch trains by default, with gaps
    between trains and injections. This is a named wrapper around
    ``fillingSchemeLHC``.

    Parameters
    ----------
    ninj : int
        Number of injections.
    nbunches : int, default 72
        Number of filled slots per train.
    ntrains : int, default 4
        Number of trains per injection.
    batchspacing : int, default 7
        Number of empty slots between trains.
    injspacing : int, default 37
        Number of empty slots after each injection.

    Returns
    -------
    list of bool
        Filling of the 3,564-slot LHC ring.
    """
    return fillingSchemeLHC(
        ninj=ninj,
        nbunches=nbunches,
        ntrains=ntrains,
        batchspacing=batchspacing,
        injspacing=injspacing,
    )


def fillingSchemeLHC_8b4e(
    ninj: int, ntrains: int = 1, batchspacing: int = 7, injspacing: int = 37
) -> list[bool]:
    """Return an LHC 8b4e filling with one or more injections.

    Each train has 56 bunches in seven groups of eight, with four empty
    slots between groups.

    Parameters
    ----------
    ninj : int
        Number of injections.
    ntrains : int, default 1
        Number of 8b4e trains per injection.
    batchspacing : int, default 7
        Number of empty slots between trains.
    injspacing : int, default 37
        Number of empty slots after each injection.

    Returns
    -------
    list of bool
        Filling of the 3,564-slot LHC ring.
    """
    # Define filling scheme: parameters
    nslots = 3564  # Defining total number of slots for SPS
    nbunches = 8 * 7  # Defining number of bunches e.g. 18, 36, 72..
    nempty = 4 * 6  # Defining number of empty slots between bunches
    # batchS = 8 # Batch spacing in 25 ns slots (200 ns)
    # injspacing = 37 # Injection spacing in 25 ns slots
    # injspacing = (nslots - ((nbunches+nempty)*ntrains+batchS*(ntrains-1))*ninj)//ninj # Uniformly distributed injection spacing

    # Defining the trains as lists of True/Falses
    bt = ([True] * 8 + [False] * 4) * 6 + [True] * 8
    st = [False] * batchspacing
    stt = [False] * injspacing
    sc = [False] * (
        nslots
        - (
            ntrains * (nbunches + nempty) * ninj
            + ((ntrains - 1) * (batchspacing) * ninj)
            + ((1) * injspacing * (ninj))
        )
    )
    an1 = (bt + st) * (ntrains - 1) + bt + stt
    an = (
        an1 * ninj + sc
    )  # This is the final true false sequence that is the beam distribution

    if len(an) > nslots:
        raise Exception(
            f"Filling scheme length: {len(an)} > available machine slots: {nslots}"
        )

    return an


# SPS filling schemes


def fillingSchemeSPS(
    nbunches: int = 72, ntrains: int = 4, batchspacing: int = 9
) -> list[bool]:
    """Return an SPS 25 ns train pattern.

    The defaults reproduce the four 72-bunch LHC25ns study pattern in
    ``SPSBeamSpectrums.ipynb``. Four 36-bunch trains model its operational
    LHC-filling beam; one bunch and one train model its AWAKE beam.

    Parameters
    ----------
    nbunches : int, default 72
        Number of consecutive filled slots in each train.
    ntrains : int, default 4
        Number of trains.
    batchspacing : int, default 9
        Number of empty 25 ns slots after each train.

    Returns
    -------
    list of bool
        Filling of the 924-slot SPS ring.
    """
    # Define filling scheme: parameters
    # ntrains = 4 # number of trains/batches
    nslots = 924  # Defining total number of slots for SPS
    # nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72..
    # batchspacing = 9 # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True] * nbunches
    st = [False] * batchspacing
    sc = [False] * (nslots - (nbunches + batchspacing) * ntrains)
    an = (bt + st) * ntrains + sc

    return an


def fillingSchemeSPS_standard(
    ntrains: int, nbunches: int = 72, batchspacing: int = 9
) -> list[bool]:
    """Return the SPS LHC25ns standard train pattern.

    By default, each train contains 72 consecutive bunches separated by
    nine empty 25 ns slots. This is a named wrapper around
    ``fillingSchemeSPS`` on the 924-slot SPS ring.

    Parameters
    ----------
    ntrains : int
        Number of trains.
    nbunches : int, default 72
        Number of filled slots per train.
    batchspacing : int, default 9
        Number of empty slots after each train.

    Returns
    -------
    list of bool
        Filling of the 924-slot SPS ring.
    """
    return fillingSchemeSPS(
        nbunches=nbunches, ntrains=ntrains, batchspacing=batchspacing
    )


def fillingSchemeSPS_BCMS(
    ntrains: int, nbunches: int = 48, batchspacing: int = 8
) -> list[bool]:
    """Return the SPS BCMS beam train pattern.

    The default train has 48 consecutive bunches followed by eight empty
    25 ns slots. This is a named wrapper around ``fillingSchemeSPS``.

    Parameters
    ----------
    ntrains : int
        Number of trains.
    nbunches : int, default 48
        Number of filled slots per train.
    batchspacing : int, default 8
        Number of empty slots after each train.

    Returns
    -------
    list of bool
        Filling of the 924-slot SPS ring.
    """
    return fillingSchemeSPS(
        nbunches=nbunches, ntrains=ntrains, batchspacing=batchspacing
    )


def fillingSchemeSPS_8b4e(ntrains: int, batchspacing: int = 8) -> list[bool]:
    """Return the SPS 8b4e beam pattern.

    Each train has six groups of eight bunches followed by four empty
    slots, then a final group of eight bunches: 56 bunches per train.

    Parameters
    ----------
    ntrains : int
        Number of trains.
    batchspacing : int, default 8
        Number of empty slots after each train.

    Returns
    -------
    list of bool
        Filling of the 924-slot SPS ring.
    """
    # Define filling scheme: parameters
    nslots = 924  # Defining total number of slots for SPS
    nbunches = 8 * 7  # Defining number of bunches e.g. 18, 36, 72..
    nempty = 4 * 6  # Defining number of empty slots between bunches
    # batchspacing = 8 # Batch spacing in 25 ns slots (200 ns)

    # Defining the trains as lists of True/Falses
    bt = ([True] * 8 + [False] * 4) * 6 + [True] * 8
    st = [False] * batchspacing
    sc = [False] * (nslots - (nbunches + nempty + batchspacing) * ntrains)
    an = (bt + st) * ntrains + sc

    return an


def fillingSchemeSPS_AWAKE(nbunches: int = 1, ntrains: int = 1) -> list[bool]:
    """Return the SPS AWAKE single-bunch pattern by default.

    The same occupancy applies before and after bunch rotation in
    ``SPSBeamSpectrums.ipynb``; the bunch shape and length belong to Beam.

    Parameters
    ----------
    nbunches : int, default 1
        Number of filled slots per train.
    ntrains : int, default 1
        Number of trains.

    Returns
    -------
    list of bool
        Filling of the 924-slot SPS ring.
    """
    # Define filling scheme: parameters
    # ntrains = 4 # number of trains/batches
    nslots = 924  # Defining total number of slots for SPS
    # nbunches = 72 # Defining a number of bunchs e.g. 18, 36, 72..
    batchspacing = 9  # Batch spacing in 25 ns slots 45/5

    # Defining the trains as lists of True/Falses
    bt = [True] * nbunches
    st = [False] * batchspacing
    sc = [False] * (nslots - (nbunches + batchspacing) * ntrains)
    an = (bt + st) * ntrains + sc

    return an


def fillingSchemeSPS_SFTPRO(
    intensities: Sequence[float], threshold: float = 1e9
) -> list[bool]:
    """Return the measured SPS SFTPRO 5 ns filling pattern.

    ``SPSBeamSpectrums.ipynb`` derives this beam from an FBCT intensity
    profile. Supply the 4,620 values on the SPS 5 ns slot grid; slots with
    intensity strictly above ``threshold`` contain a bunch.

    Parameters
    ----------
    intensities : sequence of float
        Bunch intensities in particles per slot, on the 5 ns grid.
    threshold : float, default 1e9
        Minimum intensity in particles per slot for a filled slot.

    Returns
    -------
    list of bool
        Filling of the 4,620-slot SPS SFTPRO ring.
    """
    if len(intensities) != 924 * 5:
        raise ValueError("SFTPRO intensities must contain 4620 slots")
    return [bool(intensity > threshold) for intensity in intensities]
