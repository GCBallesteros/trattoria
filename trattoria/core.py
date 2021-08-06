from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union


import numpy as np

from scipy.optimize import curve_fit

import trattoria_core


@dataclass
class G3SyncResult:
    """Results from running the synced G3 algorithm

    Attributes
    ----------
    tau1
        The tau1 values for the g3 histogram.
    tau2
        The tau2 values for the g3 histogram.
    g3
        The g3 histogram.
    """

    tau1: np.array
    tau2: np.array
    g3: np.array


@dataclass
class G3Result:
    """Results from running the G3 algorithm

    Attributes
    ----------
    tau1
        The tau1 values for the g3 histogram.
    tau2
        The tau2 values for the g3 histogram.
    g3
        The g3 histogram.
    """

    tau1: np.array
    tau2: np.array
    g3: np.array


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
class G2Parameters:
    """Parameters for the g2 algorithm."""

    channel_1: int
    channel_2: int
    correlation_window: float
    resolution: float
    record_ranges: Optional[Union[Tuple[int, int], List[Tuple[int, int]]]]


@dataclass
class ZeroFinderParameters:
    """Parameters for the zero finder algorithm.

    Attributes
    ----------
    channel_1
        The first channel of the TCSPC
    channel_2
        The second channel of the TCSPC
    correlation_window
        Correlation window length used for the fitting procedure. If not provided
        it will default to a value appropriate to the number of counts available.
    resolution
        Resolution used for the pseudo-g2 histogram. If `None` a suitable default
        is selected automatically.
    """

    channel_1: int
    channel_2: int
    correlation_window: Optional[float] = None
    resolution: Optional[float] = None


@dataclass
class LifetimeResult:
    """Results from running the lifetime algorithm

    Attributes
    ----------
    t
        The time values for the lifetime histogram.
    g2
        The lifetime histogram.
    """

    t: np.array
    g2: np.array


class TTTRFile:
    def g3(self, params: trattoria_core.G3Parameters) -> G3Result:
        """Compute the third order autocorrelation between two channels in the file.

        Arguments
        ---------
        params: G3Parameters
            A G3Parameters object has the following fields:
                - channel_1 (int): First channel used on TCSPC
                - channel_2 (int): Second channel used on TCSPC
                - channel_3 (int): Third channel used on TCSPC
                - correlation_window (float): Length in seconds of the g3 square window we
                want to consider. Values in the output histogram will be in the
                ±correlation_window range being centered at the zero delay.
                - resolution (float): Length in seconds of each bin in the g3 histogram.
                - start_record (int or None): Optionally the first record we want
                to consider when running the g3.
                - stop_record (int or None): Optionally the last record we want to
                consider when running the g3.

        Assumptions
        -----------
        1. At least three channel are being used in the TCSPC.
        2. The three channels in used are correctly identified by their integer.
        3. Failure to have similar number of counts on both channels can lead to
           asymmetries.

        Returns
        -------
        G3Result
        """
        filepath = str(self.path)
        t, hist = trattoria_core.g3(
            filepath,
            params,
        )
        hist = hist.reshape((len(t), len(t)))

        # The histogram is flipped because on the original array delays grow from
        # top to bottom and left to right but plots are usually represented with
        # axis growing towards the top-right corner.
        # We also transpose because the first index (rows) corresponds to tau1 but
        # usual representation practices put it onto the horizontal axis.
        return G3Result(tau1=np.copy(t), tau2=np.copy(t), g3=np.flipud(hist).T)

    def g3sync(self, params: trattoria_core.G3SyncParameters) -> G3SyncResult:
        """Compute the third order autocorrelation between two channels in the file.

        Arguments
        ---------
        params: G3SyncParameters
            A G3SyncParameters object has the following fields:
                - channel_sync (int): Sync channel on the the TCSPC
                - channel_1 (int): First source channel used on TCSPC
                - channel_2 (int): Second source channel used on TCSPC
                - resolution (float): Length in seconds of each bin in the g3 histogram.
                - start_record (int or None): Optionally the first record we want
                to consider when running the g3.
                - stop_record (int or None): Optionally the last record we want to
                consider when running the g3.

        Assumptions
        -----------
        1. At least three channel are being used in the TCSPC. One of them being
        a sync channel.
        2. The three channels in used are correctly identified by their integer.
        3. You are dealing with a pulsed excitation experiment

        Returns
        -------
        G3SyncResult
        """
        filepath = str(self.path)
        t, hist = trattoria_core.g3sync(
            filepath,
            params,
        )
        hist = hist.reshape((len(t), len(t)))

        # The histogram is flipped because on the original array delays grow from
        # top to bottom and left to right but plots are usually represented with
        # axis growing towards the top-right corner.
        # We also transpose because the first index (rows) corresponds to tau1 but
        # usual representation practices put it onto the horizontal axis.
        return G3SyncResult(tau1=np.copy(t), tau2=np.copy(t), g3=np.flipud(hist).T)

    def g2(self, params: G2Parameters, mode: str = "symmetric") -> G2Result:
        """Compute the second order autocorrelation between two channels in the file.

        Arguments
        ---------
        params: G2Parameters
            A G2Parameters object has the following fields:
                - channel_1 (int): First channel used on TCSPC
                - channel_2 (int): Second channel used on TCSPC
                - correlation_window (float): Length in seconds of the g2 window we
                want to consider. Values in the output histogram will be in the
                ±correlation_window range being centered at the zero delay.
                - resolution (float): Length in seconds of each bin in the g2 histogram.
                - record_range (List[Tuple[int, int]]): List of record ranges that
                should be considered in the analysis.
        mode: str
            Selects the g2 algorithm, one of "symmetric" or "asymmetric". The symmetric
            mode is more computationally expensive (x2) but gives you access to negative
            and positive time delays of channel 2 w.r.t. channel 1. In asymmetric mode
            only positive delays of channel 2 w.r.t to channel 1 are visible. In this
            case the zero delay position may not be visible in the output histogram, if
            this is the case switch the channels or use the symmetric mode.

        Assumptions
        -----------
        1. At least two channels are being used in the TCSPC.
        2. The two channels in used are correctly identified by their integer.
        3. Failure two have similar number of counts on both channels can lead to
           assymetries.

        Returns
        -------
        G2Result
        """
        filepath = str(self.path)

        # Build the trattoria_core.G2Parameters
        if params.record_ranges is not None:
            # If just have a pair of numbers instead of a list we need to wrap
            # it around a list since that is what trattoria_core expects.
            if isinstance(params.record_ranges, tuple):
                ranges = [params.record_ranges]
            else:
                ranges = params.record_ranges
        else:
            ranges = params.record_ranges

        if mode not in ["symmetric", "asymmetric"]:
            raise ValueError("Mode has to be one of symmetric or asymmetric.")

        params_core = trattoria_core.G2Parameters(
            channel_1=params.channel_1,
            channel_2=params.channel_2,
            correlation_window=params.correlation_window,
            resolution=params.resolution,
            record_ranges=ranges,
        )

        t, hist = trattoria_core.g2(
            filepath,
            params_core,
            mode=mode,
        )

        return G2Result(t=t, g2=hist)

    def timetrace(self, params: trattoria_core.TimeTraceParameters) -> TimeTraceResult:
        """Compute a time trace of the intensity as a function of time.

        The output TimeTraceResult object also contains the last record number of each
        click in a given time bin.

        Arguments
        ---------
        params: TimeTraceParameters
            A TimeTraceParameters object has the following fields:
                - resolution (float): The length of each bin of the intensity trace.
                Smaller resolutions give us a finer understanding of the intensity
                dynamics at the expense of lower counts and therefore bigger uncertainty.
                - channel (int or None): If None the intensity trace is integrated
                across channels. If provided will only include photons in said channel.

        Returns
        -------
        TimeTraceResult
        """
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

    def lifetime(self, params: trattoria_core.LifetimeParameters) -> LifetimeResult:
        """Compute the lifetime histogram from a pulsed excitation experiment.

        Arguments
        ---------
        params: LifetimeParameters
            A LifetimeParameters object has the following fields:
                - channel_sync (int): Sync channel on the TCSPC
                - channel_source (int): Channel used for the source
                - resolution (float): Length in seconds of each bin in the lifetime histogram.
                - start_record (int or None): Optionally the first record we want
                to consider when running the lifetime.
                - stop_record (int or None): Optionally the last record we want to
                consider when running the lifetime.

        Notes
        -----
        The lifetime algorithm is only supported for PTU files in T3 mode because we
        expect to have the sync rate available on the file header.

        Assumptions
        -----------
        1. We are dealing with a pulsed excitation experiment.
        2. If the sync delay line is shorter than the source line the lifetime histogram
        will appeared cut on the left.

        Returns
        -------
        LifetimeResult
        """
        filepath = str(self.path)
        t, hist = trattoria_core.lifetime(
            filepath,
            params,
        )

        return LifetimeResult(t=t, g2=hist)


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
            value, _ = self.header[field_name]
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
             t \\geq t_0 & A\\exp(\\frac{-|t-t_0|}{\\tau_r}) \\\\
           \\end{cases}
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
