# Physics Manual

This section aims to provide the user with the basic physics behind the beam induced heating computation. For a more detailed explanation, we encourage the reader to check the references provided below.

<div class="warning">
<div class="admonition-title">
Warning
</div>
This documentation is currently under development
</div>


## Introduction 
When a particle beam passes into the vacuum chamber of an accelerator component, it will electromagnetically interact with its surroundings. One of the most critical effects that occur is heating due to the passage of the beam, known as Beam-Induced Heating (BIH) [^1]. 

As a result, an accurate design stage is required to account for this effect.

This code allows the calculation of heating due to Impedance effects in a high-current accelerator ring considering:
- one single beam in the vacuum chamber
- two counter rotating beams in the same vacuum chamber

The goal of this physics guide is to introduce the definition of all the quantities that are involved in the computation of BIH and to show the mathematical formulations that are implemented in the code. 

## Wakefields and Impedances 

The longitudinal wakefunction is defined as [^2]:

$$
w_z(x_s,y_s,x_t,y_t,s) = - \frac{1}{q_s q_t} \int_0^L F_z(x_s,y_s,s,x_t,y_t,z_t)dz_t
$$

<br>

Where, $q_s$ is the charge of the source particle, in position (x_s, y_s, z_s)), $q_t$ is the charge of the test particle particle in position (x_t, y_t, z_t), $s$ is the longitudinal distance between the two particles ($s = z_s - z_t$), $L$ is the accelerator component length and, $F_z$ is the longitudinal force acting on the test particle.

The longitudinal beam-coupling impedance is defined as the Fourier transform of the longitudinal wakefunction:

$$
Z_z(\omega) = \int_{-\infty}^{\infty} w_z(s) e^{\frac{-j\omega s}{v}} \frac{ds}{v} 
$$

<br>

Where, $j$ is the imaginary unit, $\omega$ is the angular frequency conjugate variable of the time delay $\tau = \frac{s}{v}$ and, $v$ is the particle velocity.

## Single Beam and Two counter rotating beams scenarios. 
Single beam case [^3]:

$$
P_{loss} = 2 (f_0 eN_{beam})^2 \cdot \sum_{p=0}^{+\infty} |\Lambda(p\omega_0)|^2 Re[Z_z(p \omega_0)]
$$

<br>

Where, $f_0$ is the revolution frequency the machine, $e$ is the charge of the particle, $N_{beam}$ is the number of particles in the beam, $\Lambda$ is the normalized beam spectrum, $Re[Z_z]$ is the real part of the longitudinal beam-coupling impedance.

Case of two identical counter-rotating beams that are travelling in the same asymmetric vacuum chamber [^4]:

$$
P_{loss}(s) = 
(2 f_0 e N_{beam})^2 \cdot \sum_{p=0}^{+\infty} |\Lambda(p\omega_0)|^2 \cdot (Re[Z^0_z(p \omega_0)] + \\
[\Delta y_1(s) + \Delta y_2(s)]Re[Z^1_z(p \omega_0)]) \cdot (1 - cos(p \omega_0 \tau_s))
$$

<br>

Where, $s$ is the location from the interaction point at which the power loss is being computed, $Z^0_z$ and $Z^1_z$ are the longitudinal impedances of order 0 and 1, $\Delta y_1(s)$ and $\Delta y_2(s)$ are the offsets from the geometrical center of respectively Beam 1 and Beam 2, $\tau_s = 2s/c$, with $c$ the speed of light, is the relative time delay of arrival at $s$ of the two beams. 

[^1]: B. Salvant et al., “Beam induced heating”, 2012, [Online]. Available: https://cds.cern.ch/record/1975499 

[^2]: Palumbo, L., et al. "Wake Fields and Impedance". 2003. DOI.org (Datacite), https://doi.org/10.48550/ARXIV.PHYSICS/0309023.

[^3]: C. Zannini, "Electromagnetic Simulation of CERN accelerator Components and Experimental Applications", 2013. [Online]. Available: https://cds.cern.ch/record/1561199 

[^4]: C. Zannini, et al. "Power Loss Calculation in Separated and Common Beam Chambers of the LHC". Proceedings of the 5th Int. Particle Accelerator Conf., vol. IPAC2014, 2014, p. 3 pages, 1.928 MB. DOI.org (Datacite), https://doi.org/10.18429/JACOW-IPAC2014-TUPRI061. 


