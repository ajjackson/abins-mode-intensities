from itertools import groupby
from operator import itemgetter
from pathlib import Path

import dataset
from snakemake.script import snakemake

db_file = Path(snakemake.input[0])
mode_selection = snakemake.params.selection

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
        if mode_selection is not None and mode_index not in mode_selection:
            # Limit output to selected modes
            continue

        frequency = frequencies[mode_index]
        intensity = sum(map(get_intensity, group))

        print(f"{mode_index:d},{frequency:f},{intensity:f}", file=csv)
        print(f"{mode_index:4d} {frequency:8.2f} {intensity:7.3f}", file=txt)
