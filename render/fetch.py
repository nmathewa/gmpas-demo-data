"""Fetch the demo fields for 20 hourly steps of the DYAMOND2 3.75 km run
(NCAR GDEX d010129). One open per file (the open is the slow part, ~90 s),
then every field read from it. Raw float32 arrays land in OUT/<field>/."""
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import netCDF4
import numpy as np

U = "https://data.gdex.ucar.edu/d010129/dyamond2/3.75km"
OUT = Path(sys.argv[1])
FIELDS = ("olrtoa", "precipw", "refl10cm_1km", "t2m", "mslp", "uzonal_200hPa")
STAMPS = [f"2020-01-28_{h:02d}.00.00" for h in range(20)]


def get(stamp):
    todo = [f for f in FIELDS if not (OUT / f / f"{stamp}.npy").exists()]
    if not todo:
        return f"{stamp} cached"
    err = None
    for _ in range(4):
        try:
            t = time.time()
            with netCDF4.Dataset(f"{U}/diag.{stamp}.nc#mode=bytes") as nc:
                for f in todo:
                    v = np.asarray(nc[f][0, :], dtype=np.float32)
                    part = OUT / f / f"{stamp}.part.npy"
                    np.save(part, v)
                    part.rename(OUT / f / f"{stamp}.npy")
            return f"{stamp} ok {time.time() - t:.0f}s"
        except Exception as e:          # the remote read fails now and then
            err = e
            time.sleep(10)
    return f"{stamp} FAILED {err}"


if __name__ == "__main__":
    for f in FIELDS:
        (OUT / f).mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(4) as pool:  # processes: netCDF-C is not thread-safe
        for line in pool.map(get, STAMPS):
            print(line, flush=True)
    print("done", flush=True)
