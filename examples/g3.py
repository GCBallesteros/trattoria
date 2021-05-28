import os
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt

import trattoria

# G3 support is still experimental and needs to be validated against theoretical
# results.

ptu_filepath = Path("/Users/garfield/Downloads/GUI_T3_10s.ptu")
ptu = trattoria.PTUFile(ptu_filepath)
file_size = os.path.getsize(ptu_filepath) / (1024.0 ** 3)
print(ptu)

# Test G3
start_time = perf_counter()
g3_params = trattoria.G3Parameters(
    channel_1=1,
    channel_2=2,
    channel_3=3,
    correlation_window=20e-9,
    resolution=800e-12,
    start_record=None,
    stop_record=None,
)
g3_res = ptu.g3(g3_params)
end_time = perf_counter()
time_delta = end_time - start_time
print(f"G3 execution time: {time_delta:.3f} s")
print(f"  Processed {file_size/time_delta:.2f} GB/s")

t = g3_res.t * 1e9

fig = plt.figure(figsize=(8, 8))

plt.imshow(g3_res.g3, extent=(t.min(), t.max(), t.min(), t.max()))
plt.xlabel("$\\tau_1$ (ns)")
plt.ylabel("$\\tau_2$ (ns)")
plt.colorbar();

plt.show()
