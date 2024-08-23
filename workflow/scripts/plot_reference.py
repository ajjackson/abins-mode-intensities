from pathlib import Path

import matplotlib
import numpy as np
from snakemake.script import snakemake

matplotlib.use("Agg")
from matplotlib import pyplot as plt


plot_config = snakemake.params.plot
plt.style.use(plot_config["style"])


def read_spectrum(input_file: Path, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    return np.loadtxt(input_file, unpack=True, usecols=(0, 1), **kwargs)


fig, ax = plt.subplots()

x, y = read_spectrum(snakemake.input["summed"], delimiter=" ")
ax.plot(x, np.array(y) * plot_config["intensity_scale"], label="Sum over modes")

ax.plot(
    *read_spectrum(snakemake.input["fundamentals"], skiprows=2, delimiter=","),
    label="From Abins",
    linestyle=":",
)

ax.set_xlim(plot_config["x_min"], plot_config["x_max"])
ax.set_title("Fundamental-only spectrum")
ax.legend()

fig.tight_layout()
fig.savefig(snakemake.output[0])
