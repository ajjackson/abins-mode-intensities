from pathlib import Path

import mantid.simpleapi
from snakemake.script import snakemake

abins_kwargs = snakemake.params.abins_kwargs

abins_kwargs = abins_kwargs | dict(
    OutputWorkspace="fundamentals",
    SumContributions=True,
    SaveAscii=False,
    QuantumOrderEventsNumber="1",
    Autoconvolution=False,
    EnergyUnits="cm-1",
)

mantid.simpleapi.Abins(**abins_kwargs)  # Run Abins, creating new workspace
mantid.simpleapi.SaveAscii(
    InputWorkspace="fundamentals_total",
    Filename=str(Path(snakemake.output[0]).resolve()),
    AppendToFile=False,
    Separator="CSV",
)

abins_kwargs.update(
    dict(
        OutputWorkspace="multiphonon",
        QuantumOrderEventsNumber="2",
        Autoconvolution=True,
    )
)
mantid.simpleapi.Abins(**abins_kwargs)  # Run Abins, creating new workspace
mantid.simpleapi.SaveAscii(
    InputWorkspace="multiphonon_total",
    Filename=str(Path(snakemake.output[1]).resolve()),
    AppendToFile=False,
    Separator="CSV",
)
