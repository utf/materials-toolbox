"""
A script to easily access structures from the Materials Project
"""

__author__ = "Alex Ganose"
__version__ = "0.1"
__maintainer__ = "Alex Ganose"
__email__ = "alexganose@googlemail.com"
__date__ = "Feb 2, 2017"

import logging
from sys import stdout

import ase.db
from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

# ideas:
# reference


def get_structures(
    mp, query, mp_structure=False, stable=False, stable_tol=0, icsd_only=False
):
    data = ["pretty_formula", "e_above_hull", "band_gap", "nsites", "volume", "icsd_id"]
    stype = "final" if mp_structure else "initial"
    if "," in query:
        elements = query.split(",")
        entries = mp.get_entries_in_chemsys(
            elements, inc_structure=stype, property_data=data
        )
    else:
        entries = mp.get_entries(query, inc_structure=stype, property_data=data)
    if stable:
        entries = [e for e in entries if e.data["e_above_hull"] <= stable_tol]
    if icsd_only:
        entries = [e for e in entries if e.data["icsd_id"]]
    return entries


def prompt_selection(entries, save_all=False):
    from tabulate import tabulate

    headers = [
        "n",
        "Formula",
        "Spacegroup",
        "E above Hull",
        "Band Gap",
        "Nsites",
        "Volume",
        "icsd_id",
    ]
    table = []
    for i, e in enumerate(entries):
        formula = e.data["pretty_formula"]
        spg = SpacegroupAnalyzer(e.structure).get_space_group_symbol()
        e_above_hull = e.data["e_above_hull"]
        band_gap = e.data["band_gap"]
        nsites = e.data["nsites"]
        volume = e.data["volume"]
        icsd_id = e.data["icsd_id"] if e.data["icsd_id"] else "N/A"
        table.append(
            [i + 1, formula, spg, e_above_hull, band_gap, nsites, volume, icsd_id]
        )
    logging.info(tabulate(table, headers, tablefmt="orgtbl", floatfmt=".3f"))

    if save_all:
        ids = range(1, len(entries) + 1)
    else:
        id_selection = input(
            "\nWhich structures do you wish to download (type all to select all structures)?\n"
        )
        if "all" in id_selection:
            ids = range(1, len(entries) + 1)
        else:
            ids = map(int, id_selection.split())
    return ids


def save_structures(entries, ids, cif=False, conv=False, ref=False, db=None):
    """
    Args:
        ids (list): indexes of the entries to save
    """
    fmt = "cif" if cif else "poscar"
    if db:
        db = ase.db.connect(db)
    for i in ids:
        e = entries[i - 1]
        sym = SpacegroupAnalyzer(e.structure)
        if conv:
            struct = sym.get_conventional_standard_structure()
        else:
            struct = sym.get_primitive_standard_structure()

        if db is not None:
            db.write(AseAtomsAdaptor.get_atoms(struct))
        else:
            formula = e.data["pretty_formula"]
            filename = f"POSCAR_{formula}_{i}"
            struct.to(filename=filename, fmt=fmt)
        # referencing needs work
        # if ref:
        #    print e.data['icsd_id']


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="""mp-get is a convenient script that queries the
                       Materials Project database for structures""",
        epilog="""
                  Author: {}
                  Version: {}
                  Last updated: {}""".format(
            __author__, __version__, __date__
        ),
    )

    parser.add_argument(
        "query",
        help="""The material (e.g. ZnO), composition
                        (hyphen separated, e.g. Zn-O), chemical system (comma
                        separated, e.g. Zn,O) or mp-id to search for""",
    )
    parser.add_argument(
        "-s",
        "--stable",
        action="store_true",
        help="""Only include stable structures where the energy
                        above hull is 0""",
    )
    parser.add_argument(
        "--stol",
        default=0.0,
        type=float,
        help="""The tolerance in eV for stable structures
                        (defaults to 0)""",
    )
    parser.add_argument(
        "--cif",
        action="store_true",
        help="""Save structures in the cif format rather
                        than as POSCARs""",
    )
    parser.add_argument(
        "--ase",
        default=None,
        help="""The name of an ase database to save the
                        structures to (will be created if non-existant)""",
    )
    parser.add_argument(
        "--conv",
        action="store_true",
        help="""Use the conventional cell rather than the
                        primitive""",
    )
    parser.add_argument(
        "--mp_structure",
        action="store_true",
        help="""Use the Materials Project relaxed structure
                        rather than the initial ICSD structure""",
    )
    parser.add_argument(
        "--icsd_only",
        action="store_true",
        help="""Only show entries that originated from the ICSD
                        """,
    )
    parser.add_argument(
        "--ref",
        action="store_true",
        help="""Use the DOI of the reference paper as the first
                        line in the structure file""",
    )
    parser.add_argument(
        "--save_all", action="store_true", help="""Download all structures found"""
    )

    args = parser.parse_args()

    logging.basicConfig(
        filename="mp-get.log", level=logging.INFO, filemode="w", format="%(message)s"
    )
    console = logging.StreamHandler(stdout)
    logging.getLogger("").addHandler(console)

    mp = MPRester()
    entries = get_structures(
        mp,
        args.query,
        mp_structure=args.mp_structure,
        stable=args.stable,
        icsd_only=args.icsd_only,
        stable_tol=args.stol,
    )
    ids = prompt_selection(entries, save_all=args.save_all)
    save_structures(
        entries, ids, cif=args.cif, conv=args.conv, ref=args.ref, db=args.ase
    )
