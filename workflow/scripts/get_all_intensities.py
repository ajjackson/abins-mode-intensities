# import mantid algorithms, numpy
import dataset
from mantid.kernel import Atom
import mantid.simpleapi  # noqa: F401
import os
from pathlib import Path
from snakemake.script import snakemake

from abins import AbinsData, PowderCalculator
from abins.abinsalgorithm import AbinsAlgorithm
from abins.instruments import get_instrument
from abins.spowdersemiempiricalcalculator import SPowderSemiEmpiricalCalculator


def write_frequencies_to_db(frequencies, db):
    """Create a table relating mode index to frequency"""
    modes_table = db["modes"]
    for i, frequency in enumerate(frequencies):
        modes_table.insert(dict(mode_index=i, frequency=frequency))


def write_atoms_to_db(abins_data, db) -> None:
    """Create a table relating atom index to other properties"""
    atoms_table = db["atoms"]
    for i, atom in enumerate(abins_data.get_atoms_data()):
        symbol, mass = atom["symbol"], atom["mass"]
        z_number = Atom(symbol=symbol).z_number
        nucleons_number = int(round(mass))

        xs = AbinsAlgorithm.get_cross_section(
            scattering="Total", nucleons_number=nucleons_number, protons_number=z_number
        )
        atoms_table.insert(
            dict(atom_index=i, symbol=symbol, mass=mass, z=z_number, cross_section=xs)
        )


db_file = Path(snakemake.output[0])
abins_kwargs = snakemake.params.abins_kwargs
K_INDEX = snakemake.params.qpt

if db_file.exists():
    os.remove(db_file)

db = dataset.connect(f"sqlite:///{db_file.resolve()}")

# Read data from vibration/phonon calculation, write atom data to DB
abins_data = AbinsData.from_calculation_data(
    abins_kwargs["VibrationalOrPhononFile"], abins_kwargs["AbInitioProgram"]
)
write_atoms_to_db(abins_data, db)

# Set up Abins calculation objects
powder_calculator = PowderCalculator(
    filename=abins_kwargs["VibrationalOrPhononFile"],
    abins_data=abins_data,
    temperature=float(abins_kwargs["TemperatureInKelvin"]),
)
instrument = get_instrument(abins_kwargs["Instrument"], setting=abins_kwargs["Setting"])
calculator = SPowderSemiEmpiricalCalculator(
    filename=abins_kwargs["VibrationalOrPhononFile"],
    temperature=float(abins_kwargs["TemperatureInKelvin"]),
    abins_data=abins_data,
    instrument=instrument,
    quantum_order_num=1,
    autoconvolution_max=0,
)

# Compute displacement tensors and traces
powder_data = powder_calculator.get_formatted_data()
a_tensors = powder_data.get_a_tensors()[K_INDEX]
a_traces = powder_data.get_a_traces(K_INDEX)
b_tensors = powder_data.get_b_tensors()[K_INDEX]
b_traces = powder_data.get_b_traces(K_INDEX)

frequencies = powder_data.get_frequencies()[K_INDEX]
write_frequencies_to_db(frequencies, db)

# Loop over angles and atoms to get intensity contributions and write to DB
abins_table = db["abins"]
for angle in instrument.get_angles():
    q2 = instrument.calculate_q_powder(input_data=frequencies, angle=angle)

    for atom_index, atom_label in enumerate(abins_data.get_atoms_data().extract()):
        s = calculator._calculate_order_one(
            q2=q2,
            frequencies=frequencies,
            a_tensor=a_tensors[atom_index],
            a_trace=a_traces[atom_index],
            b_tensor=b_tensors[atom_index],
            b_trace=b_traces[atom_index],
        )

        dw = calculator._calculate_order_one_dw(
            q2=q2,
            frequencies=frequencies,
            a_tensor=a_tensors[atom_index],
            a_trace=a_traces[atom_index],
            b_tensor=b_tensors[atom_index],
            b_trace=b_traces[atom_index],
        )

        weights = s * dw
        for i, weight in enumerate(weights):
            abins_table.insert(
                dict(
                    angle=angle,
                    mode_index=i,
                    weight=weight,
                    atom_index=atom_index,
                    atom_label=atom_label,
                )
            )
