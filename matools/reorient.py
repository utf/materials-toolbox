"""
A script to rotate a POSCAR to align bonds on the c and x directions
Author: Alex Ganose
Version: 0.1
Date: Sept 16, 2016
"""

import ase.io
import numpy as np


def rotation_matrix(v1, v2):
    """Find rotation matrix between two vectors."""
    # Leaving this here in case you want to use a different framework without
    # an automatic rotation function
    u1 = v1 / np.linalg.norm(v1)
    u2 = v2 / np.linalg.norm(v2)
    return np.outer(u2, u1)


def angle(a, b, c):
    """Find the angle between three points."""
    ba = (a - b) / np.linalg.norm(a - b)
    bc = (c - b) / np.linalg.norm(c - b)
    angle = np.vdot(ba, bc)
    return np.arccos(angle)


def main():
    central = int(raw_input("What is the index of the central atom? ")) - 1
    along_c = int(raw_input("What is the index of the atom to align along c? ")) - 1
    along_x = int(raw_input("What is the index of the atom to align along x? ")) - 1

    atoms = ase.io.read("POSCAR")
    ang = atoms.get_angle([along_c, central, along_x]) * 180 / np.pi
    print(f"\nAngle between atoms is {ang} deg")

    # rotate so that the first bond specified is aligned along z
    # if you don't want to use ase, you can use the rotation_matrix method
    # above to find the matrix between (along_c - central) and [0, 0, 1]
    atoms.rotate(
        atoms.positions[along_c] - atoms.positions[central], [0, 0, 1], rotate_cell=True
    )

    # Because the bond angle might not always be exactly 90 deg we just rotate
    # around z until the the atom we want is above the central atom in the x
    # direction.
    # The angle to rotate is calculated by constructing two points:
    # These share the same the same x position as the along_x atom and the same
    # z position as the central atom but have different y poitions.
    central_coords = atoms.positions[central]
    along_x_coords = atoms.positions[along_x]
    p1 = np.array([along_x_coords[0], central_coords[1], central_coords[2]])
    p2 = np.array([along_x_coords[0], along_x_coords[1], central_coords[2]])
    angle = angle(p1, central_coords, p2)

    # if the rotation doesn't align the y coords of the along_x and central
    # atoms the rotation should have been in the opposite direction
    atoms.rotate("z", angle, rotate_cell=True)
    if abs(atoms.positions[along_x][1] - atoms.positions[central][1]) > 1e-4:
        angle = (2 * np.pi) - (2 * angle)
        atoms.rotate("z", angle, rotate_cell=True)

    print(f"Atom along c: {atoms.positions[along_c]}")
    print(f"Central atom: {atoms.positions[central]}")
    print(f"Atom along x: {atoms.positions[along_x]}")

    atoms.write("POSCAR_rot", vasp5=True)
    print("\nSaved rotated structure")
