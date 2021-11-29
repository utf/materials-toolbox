import argparse

"""
Script to generate the standardised primitive cell structure

Some notes:
  - The "standard" is based on: https://doi.org/10.1016/j.commatsci.2016.10.015
  - Pymatgen is used for loading files and calculating the symmetry info
  - Seekpath is used to generate the standard primitive cell as the pymatgen
    "standard" is at odds to Bradley & Cracknell and other conventions
  - Don't use spglib.find_primitive, as the find_primitive method follows a
    different convention for mC and oA as explained in the above paper
"""


def main():
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
    sym = SpacegroupAnalyzer(struct, symprec=args.tol, angle_tolerance=-1)
    data = sym.get_symmetry_dataset()

    print("Initial structure has {} atoms".format(struct.num_sites))
    print("\tSpace group number: {}".format(data["number"]))
    print("\tInternational symbol: {}".format(data["international"]))
    print("\tLattice type: {}".format(sym.get_lattice_type()))

    # first standardise the cell using the tolerance we want (seekpath has no
    # tolerance setting)
    std = spglib.refine_cell(sym._cell, symprec=args.tol)
    seek_data = seekpath.get_path(std)

    transform = seek_data["primitive_transformation_matrix"]

    # now remake the structure
    lattice = seek_data["primitive_lattice"]
    scaled_positions = seek_data["primitive_positions"]
    numbers = seek_data["primitive_types"]
    species = [sym._unique_species[i - 1] for i in numbers]
    prim = Structure(lattice, species, scaled_positions)
    prim = prim.get_sorted_structure(key=lambda x: species_order.get(x.specie.name, 0))

    Poscar(prim).write_file(f"{args.file}_prim", significant_figures=16)

    print("Final structure has {} atoms".format(prim.num_sites))
    print("Conv -> Prim transformation matrix:")
    print("\t" + str(transform).replace("\n", "\n\t"))


if __name__ == "__main__":
    main()
