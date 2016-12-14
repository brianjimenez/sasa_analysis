#!/usr/bin/env python

import timeit
import argparse
import freesasa
from freesasa import Structure
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex


def parse_command_line():
    """
    Parses command line arguments
    """
    parser = argparse.ArgumentParser(prog='scoring')
    parser.add_argument("pdb_list", help="list of PDB structures")
    script_args = parser.parse_args()
    return script_args


def calculate_sasa(pdb_file):
    # Read PDB structure
    atoms, residues, chains = parse_complex_from_file(pdb_file)
    molecule = Complex(chains, atoms, structure_file_name=pdb_file)

    # Lightdock structure to freesasa structure
    structure = Structure()
    for atom in molecule.atoms:
        structure.addAtom(atom.name, atom.residue_name, atom.residue_number, atom.chain_id, atom.x, atom.y, atom.z)

    start_time = timeit.default_timer()
    result = freesasa.calc(structure)
    elapsed = timeit.default_timer() - start_time

    return result.totalArea(), elapsed


if __name__ == "__main__":
    args = parse_command_line()

    # Set freesasa options
    freesasa.setVerbosity(freesasa.silent)

    with open(args.pdb_list) as input:
        for line in input:
            receptor, ligand = line.split()
            comp = receptor.replace('_rec.pdb', '.pdb')
            sasa_rec, time_rec = calculate_sasa(receptor)
            sasa_lig, time_lig = calculate_sasa(ligand)
            sasa_comp, time_comp = calculate_sasa(comp)
            print ','.join([str(sasa_rec), str(time_rec), str(sasa_lig), str(time_lig), str(sasa_comp), str(time_comp)])
