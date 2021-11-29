def main():
    import argparse

    from pymatgen.core.structure import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--file", type=str, default="POSCAR", help="path to input file"
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
    sym = SpacegroupAnalyzer(struct, symprec=args.tol)
    data = sym.get_symmetry_dataset()

    print("Space group number: {}".format(data["number"]))
    print("International symbol: {}".format(data["international"]))
    print(f"Lattice type: {sym.get_lattice_type()}")
