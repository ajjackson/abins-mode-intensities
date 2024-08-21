from itertools import groupby
from operator import itemgetter
from pathlib import Path

import dataset
from euphonic import ureg, Quantity, Spectrum1D
from euphonic.plot import plot_1d

import matplotlib
import numpy as np
from snakemake.script import snakemake

# Set non-interactive backend before "launching" pyplot
matplotlib.use("Agg")
import matplotlib.pyplot as plt


plot_config = snakemake.params.plot
db_file = Path(snakemake.input[0])


def tosca_broadening(x: Quantity) -> Quantity:
    """Gaussian width parameter for TOSCA as a function of energy"""
    from numpy.polynomial import Polynomial

    poly = Polynomial([2.5, 5e-3, 1e-7])
    return poly(x.to("1/cm").magnitude) * ureg("1/cm")


def plot_spectrum(x: list[float], y: list[float]) -> None:
    """Plot a binned spectrum line from scattered x/y input"""
    bins = np.linspace(plot_config["x_min"], plot_config["x_max"], 4000)
    bin_width = bins[1] - bins[0]

    hist, _ = np.histogram(x, bins=bins, weights=y, density=False)
    hist /= bin_width  # Correct intensity scaling for bins

    if snakemake.params.instrument != "TOSCA":
        raise Exception("Currently only TOSCA resolution function is available.")

    spectrum = Spectrum1D(ureg.Quantity(bins, "1/cm"), ureg.Quantity(hist, "cm"))
    spectrum = spectrum.broaden(
        x_width=tosca_broadening, shape="gauss", width_convention="std"
    )
    spectrum *= plot_config["intensity_scale"]

    fig = plot_1d(spectrum)
    fig.gca().set_title("From individual mode sums")
    fig.gca().set_xlim(bins[0], bins[-1])
    fig.savefig(snakemake.output["plot"])


# Create output directories if necessary
for output in snakemake.output:
    Path(output).parent.mkdir(parents=True, exist_ok=True)

x, y = [], []

with (
    dataset.connect(f"sqlite:///{db_file.resolve()}") as db,
    open(snakemake.output["csv"], "wt") as csv,
    open(snakemake.output["txt"], "wt") as txt,
):
    print(f"# mode_index,frequency,intensity", file=csv)
    print(f"# mode_index  frequency  intensity", file=txt)

    abins_table = db["abins"]
    atoms_table = db["atoms"]
    modes_table = db["modes"]

    cross_sections = dict(
        map(itemgetter("atom_index", "cross_section"), atoms_table.all())
    )

    frequencies = dict(map(itemgetter("mode_index", "frequency"), modes_table.all()))

    def get_intensity(data_row):
        return cross_sections[data_row["atom_index"]] * data_row["weight"]

    for mode_index, group in groupby(
        abins_table.find(order_by="mode_index"), key=itemgetter("mode_index")
    ):
        frequency = frequencies[mode_index]
        intensity = sum(map(get_intensity, group))

        print(f"{mode_index:d},{frequency:f},{intensity:f}", file=csv)
        print(f"{mode_index:4d} {frequency:8.2f} {intensity:7.3f}", file=txt)

        x.append(frequency)
        y.append(intensity)

plot_spectrum(x, y)
