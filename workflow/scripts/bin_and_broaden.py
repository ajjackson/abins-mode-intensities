from pathlib import Path

from euphonic import ureg, Quantity, Spectrum1D
import numpy as np
from numpy.polynomial import Polynomial
from snakemake.script import snakemake

bin_config = snakemake.params.bin_config


def tosca_broadening(x: Quantity) -> Quantity:
    """Gaussian width parameter for TOSCA as a function of energy"""

    poly = Polynomial([2.5, 5e-3, 1e-7])
    return poly(x.to("1/cm").magnitude) * ureg("1/cm")


def get_spectrum(x: list[float], y: list[float]) -> Spectrum1D:
    """Get a binned and broadened spectrum from scattered x/y input"""
    bins = np.linspace(bin_config["x_min"], bin_config["x_max"], bin_config["n_bins"])
    bin_width = bins[1] - bins[0]

    hist, _ = np.histogram(x, bins=bins, weights=y, density=False)
    hist /= bin_width  # Correct intensity scaling for bins

    if snakemake.params.instrument != "TOSCA":
        raise Exception("Currently only TOSCA resolution function is available.")

    spectrum = Spectrum1D(ureg.Quantity(bins, "1/cm"), ureg.Quantity(hist, "cm"))
    spectrum = spectrum.broaden(
        x_width=tosca_broadening, shape="gauss", width_convention="std"
    )
    return spectrum


def read_csv(input_file: Path) -> tuple[np.ndarray, np.ndarray]:
    return np.loadtxt(
        input_file, delimiter=",", skiprows=1, usecols=(1, 2), unpack=True
    )


spectrum = get_spectrum(*read_csv(Path(snakemake.input[0])))
spectrum.to_text_file(snakemake.output[0])
