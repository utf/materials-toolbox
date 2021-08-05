import argparse
import spglib
import seekpath

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

"""
Script to generate the standardised conventional cell structure

Some notes:
  - The "standard" is based on: https://doi.org/10.1016/j.commatsci.2016.10.015
  - Pymatgen is used for loading files and calculating the symmetry info
  - Seekpath is used to generate the standard conventional cell as the pymatgen
    "standard" is at odds to Bradley & Cracknell and other conventions
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default='POSCAR', type=str,
                        help='path to input file')
    parser.add_argument('-t', '--tol', default=1e-3, type=float,
                        help='symmetry tolerance (default 1e-3)')
    parser.add_argument('-o', '--output', default='poscar',
                        help='output file format')
    args = parser.parse_args()

    struct = Structure.from_file(args.file)
    species_order = {k.name: i for i, k in enumerate(struct.species)}
    sym = SpacegroupAnalyzer(struct, symprec=args.tol)
    data = sym.get_symmetry_dataset()

    print("Initial structure has {} atoms".format(struct.num_sites))
    print("\tSpace group number: {}".format(data['number']))
    print("\tInternational symbol: {}".format(data['international']))
    print("\tLattice type: {}".format(sym.get_lattice_type()))

    # seekpath conventional cell definition different from spglib
    std = spglib.refine_cell(sym._cell, symprec=args.tol)
    seek_data = seekpath.get_path(std)

    # now remake the structure
    lattice = seek_data['conv_lattice']
    scaled_positions = seek_data['conv_positions']
    numbers = seek_data['conv_types']
    species = [sym._unique_species[i - 1] for i in numbers]
    conv = Structure(lattice, species, scaled_positions)
    conv = conv.get_sorted_structure(key=lambda x: species_order.get(x.specie.name, 0))
    conv.to(filename="{}_conv".format(args.file), fmt=args.output)

    print("Final structure has {} atoms".format(conv.num_sites))


if __name__ == "__main__":
    main()
