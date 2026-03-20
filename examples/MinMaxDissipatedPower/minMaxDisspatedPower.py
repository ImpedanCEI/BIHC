"""
This example shows how to get the minimum and
maximum value of dissipated power when considering
a frequency uncertainty of +-20 Mhz

* date: 27/03/2023. Revised: 31/01/2024
* author: E. de la Fuente, L. Sito
"""

import sys

sys.path.append("../../")

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

import bihc

# LPC csv file name
# downloaded from LPC: https://lpc.web.cern.ch/cgi-bin/fillingSchemeTab.py
file = "25ns_2760b_2748_2494_2572_288bpi_13inj.csv"

# ------- Ploss calculation ----------

# Create beam object
bl = 1.2e-9  # bunch length [s]
Np = 2.3e11  # bunch intensity [protons/bunch]
bunchShape = "GAUSSIAN"  # bunch profile shape in time
fillMode = "FLATTOP"  # Energy
fmax = 2e9  # Maximum frequency of the beam spectrum [Hz]

beam = bihc.Beam(
    LPCfile=file,
    Np=Np,
    bunchLength=bl,
    bunchShape=bunchShape,
    machine="LHC",
    fillMode=fillMode,
    spectrum="numeric",
    fmax=fmax,
)
[f, S] = beam.spectrum

# Create Impedance object
impedance_file = "PillboxImpedance.txt"
Z = bihc.Impedance(f)
Z.getImpedanceFromCST(impedance_file)

# Get unshifted ploss
ploss, ploss_density = beam.getPloss(Z)

# ---------------- Rigid Shift power loss ------------------------------
shift = 20e6  # distance between shift steps [Hz]
shifts, power = beam.getShiftedPloss(Z, shift=shift)

print(
    f"Minimum dissipated power: P_min = {np.min(power)}, at step {shifts[np.argmin(power)]}"
)
print(
    f"Maximum dissipated power: P_max = {np.max(power)}, at step {shifts[np.argmax(power)]}"
)
print(f"Average dissipated power: P_mean = {np.mean(power)}")

# -----------   Plotting  ---------
Zmax = beam.Zmax  # Shifted impedance object
f = beam.powerSpectrum[0]  # frequency array
pow_spectrum = beam.powerSpectrum[1]  # Power spectrum Λ^2

fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(f, pow_spectrum, "b", label="Power spectrum")

axx = ax.twinx()
axx.plot(Zmax.f, Zmax.Zr / np.max(Zmax.Zr), "r", label="Shifted Impedance")
axx.plot(Z.f, Z.Zr / np.max(Z.Zr), "k", label="Impedance")
axx.set_ylabel("Normalized impedance [a.u]", color="r")
axx.set_ylim(0, 1.2)
axx.legend()

ax.text(
    0.05,
    0.95,
    f"Unshifted Power loss: {round(ploss, 2)} W",
    transform=ax.transAxes,
    color="k",
    weight="bold",
)
ax.text(
    0.05,
    0.9,
    f"Max Shifted Power loss: {round(np.max(power), 2)} W",
    transform=ax.transAxes,
    color="r",
    weight="bold",
)

ax.set_xlim(0, fmax)
ax.set_ylim(0, 1.2)
ax.set_xlabel("frequency [Hz]")
ax.set_ylabel("Power spectrum amplitude [a.u.]")
ax.set_title(
    "Maximum power loss with $f_{shift}$ of "
    + f"{round(shifts[np.argmax(power)] * shift / 1e6, 2)} MHz"
)

fig.suptitle(f'LHC power loss for Z="{impedance_file}"')
plt.show()
