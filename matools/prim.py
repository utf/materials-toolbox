import argparse
import spglib
import seekpath

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

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
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', default='POSCAR', type=str,
                        help='path to input file')
    parser.add_argument('-t', '--tol', default=1e-3, type=float,
                        help='symmetry tolerance (default 1e-3)')
    parser.add_argument('-o', '--output', default='poscar',
                        help='output file format')
    args = parser.parse_args()

    struct = Structure.from_file(args.file)
    sym = SpacegroupAnalyzer(struct, symprec=args.tol)
    data = sym.get_symmetry_dataset()

    print('Initial structure has {} atoms'.format(struct.num_sites))
    print('\tSpace group number: {}'.format(data['number']))
    print('\tInternational symbol: {}'.format(data['international']))
    print('\tLattice type: {}'.format(sym.get_lattice_type()))

    # first standardise the cell using the tolerance we want (seekpath has no
    # tolerance setting)
    std = spglib.refine_cell(sym._cell, symprec=args.tol)
    seek_data = seekpath.get_path(std)

    transform = seek_data['primitive_transformation_matrix']

    # now remake the structure
    lattice = seek_data['primitive_lattice']
    scaled_positions = seek_data['primitive_positions']
    numbers = seek_data['primitive_types']
    species = [sym._unique_species[i - 1] for i in numbers]
    prim = Structure(lattice, species, scaled_positions)
    prim.get_sorted_structure().to(filename='{}_prim'.format(args.file),
                                   fmt=args.output)

    print('Final structure has {} atoms'.format(prim.num_sites))
    print('Conv -> Prim transformation matrix:')
    print('\t' + str(transform).replace('\n', '\n\t'))
    
    # check if POSCAR atomic ordering has changed
    with open(args.file) as f: # original POSCAR file
        lines=f.readlines()
        orig_atoms_list = lines[5].split()
    
    with open('{}_prim'.format(args.file)) as f: # output prim POSCAR file
        lines=f.readlines()
        prim_atoms_list = lines[5].split()
    
    if orig_atoms_list != prim_atoms_list:
        print("Beware! POSCAR atomic ordering in POSCAR_prim does not match that of POSCAR")
        print("Make sure to update your POTCAR to match the atomic ordering of POSCAR_prim!")

if __name__ == "__main__":
    main()
