import os
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt

import trattoria

# ptu_filepath = Path("/Users/garfield/Downloads/20191205_Xminus_0p1Ve-6_CW_HBT.ptu")
ptu_filepath = Path("/Users/garfield/Downloads/GUI_T3_10s.ptu")
# ptu_filepath = Path("/Users/garfield/Downloads/GUI_T2.ptu")
ptu = trattoria.PTUFile(ptu_filepath)
file_size = os.path.getsize(ptu_filepath) / (1024.0 ** 3)
print(ptu)

# Test G2
start_time = perf_counter()
g2_params = trattoria.G2Parameters(
    channel_1=1,
    channel_2=2,
    correlation_window=10000e-12,
    resolution=60e-12,
    record_ranges=None,
)
g2_res = ptu.g2(g2_params)
end_time = perf_counter()
time_delta = end_time - start_time
print(f"G2 execution time: {time_delta:.3f} s")
print(f"  Processed {file_size/time_delta:.2f} GB/s")

plt.plot(g2_res.t * 1e9, g2_res.g2)
plt.xlabel("Delay (ns)")
plt.show()
