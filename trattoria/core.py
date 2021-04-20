from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

import numpy as np

from scipy.optimize import curve_fit

import trattoria_core


@dataclass
class G2Result:
    """Results from running the G2 algorithm

    Attributes
    ----------
    t
        The time values for the g2 histogram.
    g2
        The g2 histogram.
    """
    t: np.array
    g2: np.array


@dataclass
class TimeTraceResult:
    """Results from the time trace algorithm.

    Attributes
    ----------
    t
        The time values at which the intensity is sampled.
    tt
        The intensity time trace.
    recnum
        The last record number on each of the bins of the intensity trace.
    """
    t: np.array
    tt: np.array
    recnum: np.array


@dataclass
class ZeroFinderResult:
    """Results from the zero finder algorithm

    Attributes
    ----------
    raw_hist
        The pseudo-g2 histogram we use to derive the zero point of the histogram.
    fit_y
        The values of the fitted function used to find the zero point.
    t
        The time values corresponding to `raw_hist` and `fit_y`
    t0
        The zero point of the histogram which corresponds to the delay between
        the lines connecting to the channels of the TCSPC.
    tau1
        The tau1 parameter for the fit.
    tau2
        The tau2 parameter for the fit.
    max_intensity
        The max_intensity parameter for the fit.
    """
    raw_hist: np.array
    fit_y: np.array
    t: np.array
    t0: float
    tau1: float
    tau2: float
    max_intensity: float


@dataclass
class ZeroFinderParameters:
    channel_1: int
    channel_2: int
    correlation_window: Optional[float] = None
    resolution: Optional[float] = None


class TTTRFile:
    def g2(self, params: trattoria_core.G2Parameters):
        filepath = str(self.path)
        t, hist = trattoria_core.g2(
            filepath,
            params,
        )

        return G2Result(t=t, g2=hist)

    def timetrace(self, params: trattoria_core.TimeTraceParameters):
        filepath = str(self.path)
        timetrace, recnum_trace = trattoria_core.timetrace(
            filepath,
            params,
        )

        t = np.arange(len(timetrace)) * params.resolution

        return TimeTraceResult(t=t, tt=timetrace, recnum=recnum_trace)

    def zerofinder(self, params: ZeroFinderParameters):
        filepath = str(self.path)

        # First we need a timetrace to estimate the optimal parameters for the
        # zero finder and to seed the optimization.
        tt_params = trattoria_core.TimeTraceParameters(resolution=10, channel=None)
        tt_res = self.timetrace(tt_params)
        avg_intensity = np.mean(tt_res.tt) / tt_params.resolution

        if params.resolution is not None:
            resolution = params.resolution

        if params.correlation_window is not None:
            correlation_window = params.correlation_window

        t, zf_hist = trattoria_core.zerofinder(
            filepath,
            trattoria_core.ZeroFinderParameters(
                channel_1=params.channel_1,
                channel_2=params.channel_2,
                correlation_window=2 / avg_intensity * 4,
                resolution=(2 / avg_intensity) / 8000,
            ),
        )

        # Estimate seed parameters
        t0 = t[np.argmax(zf_hist)]
        tau = 2 / avg_intensity
        max_intensity = np.max(zf_hist)

        x_opt = fit_double_decay(t, zf_hist, tau, t0, max_intensity)

        return ZeroFinderResult(
            raw_hist=zf_hist,
            fit_y=double_decay(t, *x_opt),
            t=t,
            tau1=x_opt[0],
            tau2=x_opt[1],
            t0=x_opt[2],
            max_intensity=x_opt[3],
        )


class PTUFile(TTTRFile):
    """PTUFile contain the metadata for a PicoQuant PTU file and the streaming algorithms."""

    def __init__(self, filepath: Path):
        """Initialize a PTUFile.

        The initialized object contains the file metadata within  the `self.header` dict.
        `self.header` each value in self.header is a tuple that contains the value and
        a string describing the value type as read from the file header.
        """
        self.path = filepath
        self.header = trattoria_core.read_ptu_header(str(self.path))

    def __str__(self) -> str:
        """Pretty print the metadata dictionary to the screen."""
        fields: List[str] = []
        for field_name in self.header:
            value, ftype = self.header[field_name]
            fields.append(f"{field_name:<35} {value}")

        nl = "\n"
        return f"{self.path}{nl}{nl.join(fields)}"

    def __repr__(self) -> str:
        return (repr(self.ptu), repr(self.header))


def double_decay(t, tau1, tau2, t0, max_value):
    """Asymmetricly decaying exponential.

    ```
    f(t) = \\begin{cases} 
             t<t_0 & A\exp(\\frac{-|t-t_0|}{\\tau_l}) \\\\
             t \geq t_0 & A\exp(\\frac{-|t-t_0|}{\\tau_r}) \\\\
           \end{cases}
    ```

    Parameters
    ----------
    t
        The times at which we want to compute the function.
    tau1
        The decay constant for the left part of the function.
    tau2
        The decay constant for the right part of the function.
    t0
        Center position of the function
    max_value
        Value at t0

    Returns
    -------
    Array with the values of the function.
    """
    tau = np.where(t < t0, tau1, tau2)

    decay = max_value * np.exp(-np.abs((t - t0) / tau))

    return decay


def fit_double_decay(t, y, tau, t0, max_value):
    x_opt = curve_fit(double_decay, t, y, p0=[tau, tau, t0, max_value])[0]

    return x_opt
