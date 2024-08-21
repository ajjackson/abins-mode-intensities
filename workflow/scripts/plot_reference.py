from pathlib import Path

import matplotlib
import numpy as np
from snakemake.script import snakemake

matplotlib.use("Agg")
from matplotlib import pyplot as plt


plot_config = snakemake.params.plot


def read_spectrum(input_file: Path) -> tuple[np.ndarray, np.ndarray]:
    return np.loadtxt(
        input_file, delimiter=",", skiprows=2, usecols=(0, 1), unpack=True
    )


fig, ax = plt.subplots()

ax.plot(*read_spectrum(snakemake.input["multiphonon"]), label="Multi-phonon")
ax.plot(*read_spectrum(snakemake.input["fundamentals"]), label="Fundamentals")

ax.set_xlim(plot_config["x_min"], plot_config["x_max"])
ax.set_title("Abins total spectrum")
ax.legend()

fig.tight_layout()
fig.savefig(snakemake.output[0])
