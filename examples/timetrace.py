import os
from pathlib import Path
from time import perf_counter

import matplotlib.pyplot as plt

import trattoria

# ptu_filepath = Path("/Users/garfield/Downloads/20191205_Xminus_0p1Ve-6_CW_HBT.ptu")
ptu_filepath = Path("/Users/garfield/Downloads/GUI_T3_10s.ptu")
ptu = trattoria.PTUFile(ptu_filepath)
file_size = os.path.getsize(ptu_filepath) / (1024.0 ** 3)
print(ptu)

start_time = perf_counter()
timetrace_params = trattoria.TimeTraceParameters(
    resolution=1.0,
    channel=4,
)
tt_res = ptu.timetrace(timetrace_params)
end_time = perf_counter()
time_delta = end_time - start_time

print(f"Timetrace execution time: {time_delta:.3f} s")
print(f"  Processed {file_size/time_delta:.2f} GB/s")

plt.plot(tt_res.t, tt_res.tt)
plt.xlabel("Time (s)")
plt.show()
