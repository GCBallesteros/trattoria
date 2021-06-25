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

# Test G2 With manual postselection
start_time = perf_counter()


# We can manually specify what records to consider.
# Note the explicit casting to integers. G2Parameters expects a pair of actual
# integers, hence the int(1e7) and the use of tuples. For example the following
# is not valid:
# [[0, int(1e7)]] <- invalid. This a list of lists not a list of tuples
# [(0, 1e7)] <- invalid. This a list of tuples but one of the values is a float!
manual_post_selection_ranges = [
    (0, int(1e7)),
    (int(3e7), int(4e7)),
]

# you can also derive the postselection intervals based on an intensity threshold
timetrace_params = trattoria.TimeTraceParameters(
    resolution=0.1,
    channel=0,
)
tt_res = ptu.timetrace(timetrace_params)

auto_post_selection_ranges = trattoria.construct_postselect_vector(
    tt_res,  # Result from a timetrace experiment
    1000,  # Counts per time resolution bin in the trace
    True,  # True = select above the threshold, False = select below
)

g2_params = trattoria.G2Parameters(
    channel_1=1,
    channel_2=2,
    correlation_window=10000e-12,
    resolution=60e-12,
    record_ranges=auto_post_selection_ranges,
)
g2_res = ptu.g2(g2_params)
end_time = perf_counter()
time_delta = end_time - start_time
print(f"G2 execution time: {time_delta:.3f} s")
print(f"  Processed {file_size/time_delta:.2f} GB/s")

plt.plot(g2_res.t * 1e9, g2_res.g2)
plt.xlabel("Delay (ns)")
plt.show()
