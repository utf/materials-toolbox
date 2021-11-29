"""
Script to generate the standardised conventional cell structure

Some notes:
  - The "standard" is based on: https://doi.org/10.1016/j.commatsci.2016.10.015
  - Pymatgen is used for loading files and calculating the symmetry info
  - Seekpath is used to generate the standard conventional cell as the pymatgen
    "standard" is at odds to Bradley & Cracknell and other conventions
"""


def main():
    import argparse

    import seekpath
    import spglib
    from pymatgen.core.structure import Structure
    from pymatgen.io.vasp.inputs import Poscar
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--file", default="POSCAR", type=str, help="path to input file"
    )
    parser.add_argument(
        "-t",
        "--tol",
        default=1e-3,
        type=float,
        help="symmetry tolerance (default 1e-3)",
    )
    args = parser.parse_args()

    struct = Structure.from_file(args.file)
    species_order = {k.name: i for i, k in enumerate(struct.species)}
    sym = SpacegroupAnalyzer(struct, symprec=args.tol)
    data = sym.get_symmetry_dataset()

    print(f"Initial structure has {struct.num_sites} atoms")
    print("\tSpace group number: {}".format(data["number"]))
    print("\tInternational symbol: {}".format(data["international"]))
    print(f"\tLattice type: {sym.get_lattice_type()}")

    # seekpath conventional cell definition different from spglib
    std = spglib.refine_cell(sym._cell, symprec=args.tol)
    seek_data = seekpath.get_path(std)

    # now remake the structure
    lattice = seek_data["conv_lattice"]
    scaled_positions = seek_data["conv_positions"]
    numbers = seek_data["conv_types"]
    species = [sym._unique_species[i - 1] for i in numbers]
    conv = Structure(lattice, species, scaled_positions)
    conv = conv.get_sorted_structure(key=lambda x: species_order.get(x.specie.name, 0))

    Poscar(conv).write_file(f"{args.file}_conv", significant_figures=16)

    print(f"Final structure has {conv.num_sites} atoms")
