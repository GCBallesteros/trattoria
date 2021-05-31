import os
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt

import numpy as np

import trattoria

ptu_filepath = Path("/Users/garfield/Downloads/GUI_T3_10s.ptu")
# ptu_filepath = Path("/Users/garfield/Downloads/28102020_QD_TemporalPurity_Test1_match_arrival_004.ptu")
ptu = trattoria.PTUFile(ptu_filepath)
file_size = os.path.getsize(ptu_filepath) / (1024.0 ** 3)
print(ptu)  # Prints the header

# Test G3
start_time = perf_counter()
g3_params = trattoria.G3SyncParameters(
    channel_sync=0,
    channel_1=1,
    channel_2=2,
    resolution=80e-12,
    start_record=None,
    stop_record=None,
)
g3_res = ptu.g3sync(g3_params)
end_time = perf_counter()
time_delta = end_time - start_time
print(f"G3 execution time: {time_delta:.3f} s")
print(f"  Processed {file_size/time_delta:.2f} GB/s")

tau1 = g3_res.tau1 * 1e9
tau2 = g3_res.tau2 * 1e9

fig = plt.figure(figsize=(8, 8))

plt.imshow(np.log10(g3_res.g3), extent=(tau1.min(), tau1.max(), tau2.min(), tau2.max()))
plt.xlabel("$\\tau_1$ (ns)")
plt.ylabel("$\\tau_2$ (ns)")
plt.colorbar();

plt.show()
