from pathlib import Path

import matplotlib
import numpy as np
from snakemake.script import snakemake

# Set non-interactive backend before "launching" pyplot
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots()

plot_config = snakemake.params.plot
plt.style.use(plot_config["style"])


def read_spectrum(input_file: Path, **kwargs) -> tuple[np.ndarray, np.ndarray]:
    return np.loadtxt(input_file, unpack=True, usecols=(0, 1), **kwargs)


fig, ax = plt.subplots()

# Plot multiphonon spectrum from Abins
ax.plot(
    *read_spectrum(snakemake.input["multiphonon"], delimiter=",", skiprows=2),
    label="Fundamentals + multiphonon"
)

# Plot fundamental spectrum summed by Abins
ax.plot(
    *read_spectrum(snakemake.input["fundamentals"], delimiter=",", skiprows=2),
    label="Fundamentals"
)

# Plot contribution from selected modes
x, y = read_spectrum(snakemake.input["highlights"], delimiter=" ")
ax.plot(x, y * plot_config["intensity_scale"], label="Selected modes")


ax.set_title("Simulated INS spectrum")
ax.legend()
ax.set_xlabel("Frequency / cm$^{-1}$")
ax.set_xlim(plot_config["x_min"], plot_config["x_max"])
ax.set_ylabel("Intensity")
fig.tight_layout()
fig.savefig(snakemake.output["plot"])
