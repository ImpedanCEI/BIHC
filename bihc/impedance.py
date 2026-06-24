# copyright ################################# #
# This file is part of the BIHC Package.      #
# Copyright (c) CERN, 2024.                   #
# ########################################### #

"""
Impedance module to manage impedance object creation from a txt file exported
from CST database or thanks to (i) the resisitive wall impedance model (for a
circular pipe) or (ii) broad-band resonator model.

The created impedance consists in a list of numpy.ndarray: the first element is
the frequency vector, the second is the impedance.

* date: 12/12/2022
* author: Francesco Giordano, Elena de la Fuente, Leonardo Sito
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import NDArray

from bihc.plot import Plot


class Impedance(Plot):
    """Define a Beam Coupling Impedance curve.

    It is possible to define an impedance curve using: (i) the broad-band
    resonator model, (ii) the resistive wall impedance of a circular beam pipe,
    and (iii) importing an external impedance text file (like the one that can
    be exported from CST).

    Mixin class for the Beam class.

    Parameters
    ----------
    f : numpy.ndarray, default np.linspace(0.1,2e9,int(1e5))
        Frequency array of the impedance curve.

        This is the numpy array containing the frequency point of the impedance
        curve. by default it is a linearly spaced vector going from 0.1 Hz to 2
        GHz with 10'000 points.

    CST_file : str, deafult 0
        File name with extension of the txt file containing the impedance curve.

    Attributes
    ----------
    f : numpy.ndarray, default np.linspace(0.1,2e9,int(1e5))
        Frequency array of the impedance curve [Hz].
    Zr : numpy.ndarray
        Real part of the impedance in the speciefied frequency points [Ohm].
    Zi : numpy.ndarray
        Imaginary part of the impedance in the speciefied frequency points [Ohm].
    """

    def __init__(
        self,
        f: NDArray[np.floating[Any]] = np.linspace(0.1, 2e9, int(1e5)),
        Z: NDArray[np.complexfloating[Any, Any]] | None = None,
        CST_file: str | None = None,
    ) -> None:
        self.f = f
        self.Zr = np.zeros(len(f))
        self.Zi = np.zeros(len(f))

        if Z is not None:
            self.Zr = np.real(Z)
            self.Zi = np.imag(Z)
            if len(self.f) != len(self.Zr):
                print("[!] frequency array and impedance data have different lengths")

        self.isResonatorImpedance = False
        self.isRWImpedance = False

        if CST_file is not None:
            self.getImpedanceFromCST(CST_file)

    def __add__(self, Zn: Impedance) -> Impedance:
        Z = Impedance(self.f)

        if not np.array_equal(self.f, Zn.f):
            Z.Zr = np.interp(self.f, Zn.f, Zn.Zr)
            Z.Zi = np.interp(self.f, Zn.f, Zn.Zi)
        else:
            Z.Zr = Zn.Zr
            Z.Zi = Zn.Zi

        Z.Zr += self.Zr
        Z.Zi += self.Zi

        return Z

    def getResonatorImpedance(
        self, Rs: float, Qr: float, fr: float
    ) -> list[NDArray[np.floating[Any]]]:
        """Creating the impedance curve from the broad-band resonator model.

        This methods creates an impedance curve using the broad-band resonator
        model. It requires a shunt impedance value, a quality factor value, and
        the resonant frequency value.

        Parameters
        ----------
        Rs : float
            Shunt impedance in the broad-band resonator model [Ohm].
        Qr : float
            Quality factor in the broad-band resonator model.
        fr : float
            Resonant frequency in the broad-band resonator model [Hz].

        Returns
        -------
        [f, Z] : numpy.ndarray list
            Impedance curve. Returns the list of numpy arrays [frequency,
            complex impedance].
        """

        self.fr = fr
        self.Rs = Rs
        self.Qr = Qr
        self.isResonatorImpedance = True

        f = self.f
        mask1 = f == 0
        f[mask1] = 1e-5
        Z = Rs / (1 + 1j * Qr * (f / fr - fr / f))

        self.Zr = np.real(Z)
        self.Zi = np.imag(Z)
        self.Z = Z

        return [self.f, self.Zr + 1j * self.Zi]

    def getRWImpedance(
        self, L: float, b: float, sigma: float
    ) -> list[NDArray[np.floating[Any]]]:
        """Creating the impedance curve from the resistive wall impedance model.

        This methods creates an impedance curve using the resistive wall
        impedance model for a resistive pipe of circular section with.
        It considers the thick wall regime.

        Parameters
        ----------
        L : float
            Length of the pipe in the resistive wall impedance model [m].
        b : float
            Pipe radius in the resisitve wall impedance model [m].
        sigma : float
            Electrical conductivity of the pipe in the resistive wall
            impedance model [Siemens/m].

        Returns
        -------
        [f, Z] : numpy.ndarray list
            Impedance curve. Returns the list of numpy arrays [frequency,
            complex impedance].
        """
        Z0 = 376.73
        c = 3e8
        f = self.f

        # !n.b. problem with the signum
        Z = (
            (1 + 1j * np.sign(f))
            * L
            / (2 * np.pi * b)
            * np.sqrt(Z0 * np.abs(2 * np.pi * f) / (2 * c * sigma))
        )

        self.Zr = np.real(Z)
        self.Zi = np.imag(Z)
        self.Z = Z

        self.b = b
        self.sigma = sigma
        self.L = L
        self.isRWImpedance = True

        return [self.f, self.Zr + 1j * self.Zi]

    def getImpedanceFromCST(
        self, path: str, unit: str = "GHz", skip_header: int = 2, skip_footer: int = 0
    ) -> list[NDArray[np.floating[Any]]]:
        """Creating the impedance curve from CST file.

        This methods creates an impedance curve using a txt file, usually
        exported from CST.

        Parameters
        ----------
        path : str
            File name with extension.
        unit : str, default "GHz"
            Units of the frquency array from the simulation. Either GHz or MHz.
        skip_header : int, default 2
            Line index to skip to at the beginning of the txt file.
        skip_footer : int, default 0
            Line index from which to skip to at the end of the txt file.

        Returns
        -------
        [f, Z] : numpy.ndarray list
            Impedance curve. Returns the list of numpy arrays [frequency,
            complex impedance].
        """
        data = np.genfromtxt(path, skip_header=skip_header, skip_footer=skip_footer)

        if unit == "GHz":
            self.f = data[:, 0] * 1e9
        elif unit == "MHz":
            self.f = data[:, 0] * 1e6

        self.Zr = data[:, 1]
        try:
            self.Zi = data[:, 2]
        except Exception:
            print("[!] No imaginary part found in the CST file. Setting it to zero.")
            self.Zi = np.zeros_like(self.Zr)

        self.Z = self.Zr + 1j * self.Zi

        return [self.f, self.Z]

    def getImpedanceFromFunc(
        self,
        func: Callable[
            [NDArray[np.floating[Any]]], NDArray[np.complexfloating[Any, Any]]
        ],
        f: NDArray[np.floating[Any]] | None = None,
    ) -> list[NDArray[np.floating[Any]]]:
        """Gets impedance from a python function
        in the form of Z = func(f)

        Ensures compatibility to Xwakes wake objects
        https://github.com/xsuite/xwakes/

        Example:
        -------
        >>> import xwakes as xw
        >>> wake = xw.WakeResonator(r=1e8, q=10, f_r=1e9,
                                    kind='longitudinal')
        >>> import bihc
        >>> Z = bihc.Impedance(f=np.linspace(0,1e9,10000))
        >>> Z.getImpedanceFromFunc(wake.components[0].impedance)
        """

        if f is not None:
            self.f = f

        self.Z = func(f)
        self.Zr = np.real(self.Z)
        self.Zi = np.imag(self.Z)

        return [self.f, self.Z]

    def getImpedanceFromPandas(self, path: str, unit: str = "GHz") -> Any:
        import pandas as pd

        data = pd.read_csv(path)
        data.columns = ["f", "Zr"]

        if unit == "GHz":
            data["f"] = data["f"] * 1e9
        elif unit == "MHz":
            data["f"] = data["f"] * 1e6

        self.f = data["f"]
        self.Zr = data["Zr"]

        return data

    def copy(self) -> Impedance:
        obj = type(self).__new__(self.__class__)
        obj.__dict__.update(self.__dict__)
        return obj

    def getFrequencyRegions(
        self,
        vlines: NDArray[np.floating[Any]] | None = None,
        figsize: list[float] = [12, 6],
        dpi: int = 200,
    ) -> NDArray[np.floating[Any]]:
        """Interactively select frequency points of the impedance curve.

        The picked frequencies are stored in ``self.freqregions`` and returned.
        Behaviour depends on the active Matplotlib backend:

        * **Interactive / script** (e.g. ``Qt5Agg``, ``TkAgg``):
          Spacebar --> pick | Delete --> unpick | Enter --> Finish

        * **Jupyter notebook with** ``%matplotlib widget`` **(ipympl)**:
          Left click --> pick | Right click --> undo last pick.
          A *Done* button is displayed below the figure; click it to finish.

          .. note::
             The function returns immediately with an empty array.
             Results are stored in ``self.freqregions`` when *Done* is pressed.

        * **Jupyter notebook with** ``%matplotlib inline``:
          A static plot is shown and you are prompted to type the frequencies
          (in GHz) as a space- or comma-separated list.

        Parameters
        ----------
        vlines : numpy.ndarray, optional
            Frequencies [Hz] at which to draw dashed reference lines.
        figsize : list, default [12, 6]
            Figure size ``[width, height]`` in inches.
        dpi : int, default 200
            Figure resolution.

        Returns
        -------
        freqregions : numpy.ndarray
            Selected frequencies [Hz].  In ``%matplotlib widget`` mode the
            array is empty on return and is populated once *Done* is clicked.
        """
        import matplotlib
        import matplotlib.pyplot as plt

        self.freqregions = np.array([])
        backend = matplotlib.get_backend().lower()
        is_inline = "inline" in backend
        is_widget = "ipympl" in backend or "nbagg" in backend

        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        ax.plot(self.f * 1e-9, self.Zr, color="k")

        if vlines is not None:
            ax.vlines(
                vlines * 1e-9,
                ymin=np.min(self.Zr),
                ymax=np.max(self.Zr),
                color="k",
                linestyle="dashed",
                alpha=0.4,
            )

        ax.set_ylabel(r"Impedance [$\Omega$]")
        ax.set_xlabel("Frequency [GHz]")
        fig.tight_layout()

        # ---------------------------------------------------------------- inline
        if is_inline:
            ax.set_title("Static view — enter frequencies below", color="red")
            plt.show()
            response = input(
                "Enter frequencies in GHz separated by spaces or commas: "
            ).strip()
            if response:
                vals = response.replace(",", " ").split()
                self.freqregions = np.array([float(v) for v in vals]) * 1e9

            plt.close(fig)
            print(f"Selected frequency regions: {self.freqregions} Hz")
            return self.freqregions

        # ---------------------------------------------------- notebook widget
        if is_widget:
            ax.set_title(
                "Left click to pick | Right click to undo | Click Done when finished",
                color="red",
            )
            picked = []
            vline_artists = []

            def on_click(event):
                if event.inaxes != ax or event.xdata is None:
                    return
                if event.button == 1:  # left click → add
                    picked.append(event.xdata)
                    vl = ax.axvline(event.xdata, color="r", linestyle="--", alpha=0.7)
                    vline_artists.append(vl)
                    fig.canvas.draw_idle()
                elif event.button == 3:  # right click → undo last
                    if picked:
                        picked.pop()
                        vline_artists.pop().remove()
                        fig.canvas.draw_idle()

            fig.canvas.mpl_connect("button_press_event", on_click)

            try:
                import ipywidgets as widgets
                from IPython.display import display as ipy_display

                status_label = widgets.Label(value="Selected: none")
                done_btn = widgets.Button(
                    description="Done", button_style="success", icon="check"
                )
                out = widgets.Output()

                def _update_status():
                    if picked:
                        entries = ", ".join(f"{p:.4f} GHz" for p in picked)
                        status_label.value = f"Selected: {entries}"
                    else:
                        status_label.value = "Selected: none"

                # Patch on_click to update the live label after each pick/undo
                _orig_on_click = on_click

                def on_click(event):  # noqa: F811
                    _orig_on_click(event)
                    _update_status()

                fig.canvas.mpl_connect("button_press_event", on_click)

                def on_done(b):
                    self.freqregions = np.array(picked) * 1e9
                    done_btn.disabled = True
                    status_label.value = (
                        f"Done — {len(self.freqregions)} region(s) stored in "
                        "self.freqregions"
                    )
                    with out:
                        print(
                            f"Selected frequency regions: {self.freqregions * 1e-9} GHz"
                        )
                    plt.close(fig)

                done_btn.on_click(on_done)
                plt.show()
                ipy_display(widgets.VBox([status_label, done_btn, out]))
                print(
                    "Click on the plot to pick frequencies. "
                    "Access self.freqregions after clicking Done."
                )
            except ImportError:
                plt.show()
                print(
                    "Click on the plot to pick frequencies. "
                    "Close the figure when done — "
                    "results will be stored in self.freqregions."
                )

                def on_close(event):
                    self.freqregions = np.array(picked) * 1e9
                    print(f"Selected frequency regions: {self.freqregions * 1e-9} GHz")

                fig.canvas.mpl_connect("close_event", on_close)

            # Kernel returns immediately; self.freqregions is populated on Done
            self.freqregions = np.array([])
            return self.freqregions

        # ------------------------------------------------------- desktop / .py
        ax.set_title(
            "Spacebar to pick | Delete to unpick | Enter to finish",
            color="red",
        )
        points = plt.ginput(
            n=0, timeout=0, mouse_add=None, mouse_pop=None, mouse_stop=None
        )
        points = np.asarray(points, dtype=float)
        self.freqregions = points[:, 0] * 1e9 if points.size else np.array([])
        plt.close(fig)
        print(f"Selected frequency regions: {self.freqregions} Hz")
        return self.freqregions
