import sys
import argparse

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dim', type=str, nargs='+',
                        help="The supercell dimensions, either as 3 numbers, "
                              "a full 3x3 scaling matrix, or a single scaling "
                              "factor")
    parser.add_argument('-f', '--file', type=str, default='POSCAR',
                        help='path to input file')
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
    nsites = struct.num_sites
    struct.make_supercell(dim)
    struct.to(filename="{}_super".format(args.file), fmt="poscar")

    print("Initial structure has {} atoms".format(nsites))
    print("Final structure has {} atoms".format(struct.num_sites))

if __name__ == "__main__":
    main()
