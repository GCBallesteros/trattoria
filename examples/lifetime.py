import os
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt

import trattoria

ptu_filepath = Path("/Users/garfield/Downloads/GUI_T3_10s.ptu")
ptu = trattoria.PTUFile(ptu_filepath)
file_size = os.path.getsize(ptu_filepath) / (1024.0 ** 3)
print(ptu)

# Test Lifetime
start_time = perf_counter()
lifetime_params = trattoria.LifetimeParameters(
    channel_sync=0,
    channel_source=2,
    resolution=60e-12,
    start_record=None,
    stop_record=None,
)
lifetime_res = ptu.lifetime(lifetime_params)
end_time = perf_counter()
time_delta = end_time - start_time
print(f"Lifetime execution time: {time_delta:.3f} s")
print(f"  Processed {file_size/time_delta:.2f} GB/s")

plt.plot(lifetime_res.t * 1e9, lifetime_res.g2)
plt.xlabel("Delay (ns)")
plt.show()
