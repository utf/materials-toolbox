import sys
import argparse

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', type=str, default='POSCAR',
                        help='path to input file')
    parser.add_argument('-t', '--tol', default=1e-3, type=float,
                        help='symmetry tolerance (default 1e-3)')
    args = parser.parse_args()

    struct = Structure.from_file(args.file)
    sym = SpacegroupAnalyzer(struct, symprec=args.tol)
    data = sym.get_symmetry_dataset()

    print("Space group number: {}".format(data['number']))
    print("International symbol: {}".format(data['international']))
    print("Lattice type: {}".format(sym.get_lattice_type()))


if __name__ == "__main__":
    main()
