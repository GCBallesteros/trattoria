# 🍕 Trattoria 🍕
Trattoria delivers you the fastest streaming algorithms to analyze your TTTR data. We
currenlty support the following algorithms:
- __Second order autocorrelations__: Calculate the autocorrelation between two channels of
  your TCSPC.
- __Third Order autocorrelations__: Calculate the coincidences between 3 channels. A sync
version is provided were it uses the fact that the sync channel is periodic and known.
- __Intensity time trace__: Calculate the intensity on each (or all) channels versus time.
- __Zero finder__: Given two uncorrelated channels (e.g. a laser behind a 50/50 splitter)
  compute the delay between the input channels.
- __Lifetime__: Compute the lifetime histogram from a pulsed excitation experiment.

## Supported file formats
Currently Trattoria can only read PTU files from PicoQuant. If you want support
for more or want to help providing it please put a ticket on the tttr-toolbox
project.

## Installing
```
pip install trattoria
```

## Examples
The entry point to Trattoria is the PTUFile class. This class has methods that
give us access to the algorithms. Each of the algorithms takes as input a
parameter object and returns a results object. For a complete list of the
functionality see the `examples` folder.

```python
from pathlib import Path

import trattoria

import matplotlib.pyplot as plt

ptu_filepath = Path("/path/to/some.ptu")
ptu = trattoria.PTUFile(ptu_filepath)

timetrace_params = trattoria.TimeTraceParameters(
    resolution=10.0,
    channel=None,
)
tt_res = ptu.timetrace(timetrace_params)

plt.plot(tt_res.t, tt_res.tt / timetrace_params.resolution)
plt.xlabel("Time (s)")
plt.ylabel("Intensity (Hz)")
plt.show()
```

The examples folders contains examples of all the functionality available in
Trattoria.  For more details check the docstrings in `core.py`.

## Design
Trattoria is just a very thin wrapper around the
[trattoria-core](https://github.com/GCBallesteros/trattoria-core) library which
itselfs provides a lower level interface to the the
[tttr-toolbox](https://github.com/GCBallesteros/tttr-toolbox/tree/master/tttr-toolbox)
library. A Rust project that provides the compiled components that allows us to
go fast.

## Changelog
### 0.3.5
- Bug fix. The last 1024\*16 where being ignored for performance reasons. This has
  has been fixed upstream in `tttr-toolbox` and this version of Trattoria uses the
  upgraded version of `trattoria-core`.
- `trattoria-core` dropped support for Python 3.6 and 3.7 and therefore Trattoria too.

### 0.3.4
- The g2 algorithm now supports a mode flag. With "symmetric" we use the
  prefered version of the algorithm that returns negative and positive delays.
  "asymmetric" returns only positive delays but is faster. Default is
  "symmetric".

### 0.3.3
- The underlying TTTR Toolbox and Trattoria Core were refactored to support
  multiple custom ranges or records at once. `start_range` and `stop_range`
  have disappeared in favor of `record_ranges`. It takes a list of tuples of
  integers or `None`.


## Citing

