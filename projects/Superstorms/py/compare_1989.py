import sys

import pandas as pd

sys.path.append("py/")
import datetime as dt

from plots import TimeSeriesPlot

data = pd.read_csv(
    "../ScaleHistoricData/March1989/Voltage/TAT8Volt-rescale.csv", parse_dates=["Time"]
)
data.rename(columns={"Voltage": "Vt(v)"}, inplace=True)
data.sort_values("Time", inplace=True)
data.set_index("Time", inplace=True)
model = pd.read_csv("../../1989Storm.csv", parse_dates=["Time"])
model.set_index("Time", inplace=True)
ts = TimeSeriesPlot(
    [dt.datetime(1989, 3, 13, 12), dt.datetime(1989, 3, 14, 12)],
    "",
    num_subplots=1,
    text_size=15,
)
model["Vt(v)"] = -model["Vt(v)"]  # Why need to check with David.
# print(data.tail(), data.head())
ax = ts.add_voltage(model, xlabel="")
ts.add_voltage(data, color="r", xlabel="Time since 12 UT, 12 March 1989 [UT]", ax=ax)
ts.save("figures/compare.png")
ts.close()
