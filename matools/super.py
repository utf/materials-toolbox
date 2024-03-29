def main():
    import argparse
    import sys

    import numpy as np
    from pymatgen.core.structure import Structure
    from pymatgen.io.vasp.inputs import Poscar

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dim",
        type=str,
        nargs="+",
        help="The supercell dimensions, either as 3 numbers, "
        "a full 3x3 scaling matrix, or a single scaling "
        "factor",
    )
    parser.add_argument(
        "-f", "--file", type=str, default="POSCAR", help="path to input file"
    )
    args = parser.parse_args()

    int_dim = list(map(int, args.dim))
    if len(int_dim) == 1:
        dim = int_dim[0]
    elif len(int_dim) == 3:
        dim = int_dim
    elif len(int_dim) == 9:
        dim = np.array(int_dim).reshape(3, 3)
    else:
        print("Cannot parse supercell dimensions, try `super -h` for help")
        sys.exit()

    struct = Structure.from_file(args.file)
    species_order = {k.name: i for i, k in enumerate(struct.species)}
    nsites = struct.num_sites
    struct.make_supercell(dim)
    struct = struct.get_sorted_structure(
        key=lambda x: species_order.get(x.specie.name, 0)
    )
    struct.to(filename=f"{args.file}_super", fmt="poscar")
    Poscar(struct).write_file(f"{args.file}_super", significant_figures=16)

    print(f"Initial structure has {nsites} atoms")
    print(f"Final structure has {struct.num_sites} atoms")
