#!/usr/bin/env python

import timeit
import argparse
import access
import numpy as np
from lightdock.pdbutil.PDBIO import parse_complex_from_file
from lightdock.structure.complex import Complex
from libpydock.parameters.DesolvationParameters import DesolvationParameters


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

    # Load parameters
    parameters = DesolvationParameters()

    x = []
    y = []
    z = []
    atom_names = []
    atom_radius = []
    for atom in molecule.atoms:
        x.append(atom.x)
        y.append(atom.y)
        z.append(atom.z)
        atom_names.append("%-4s" % atom.name)
        if atom.residue_name == 'CYX':
            atom.residue_name = 'CYS'
        atom_radius.append(parameters.radius_per_atom[atom.residue_name + "-" + atom.name])
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    num_residues = len(molecule.residues)
    num_atoms = len(molecule.atoms)
    residue_ids = [residue.number for residue in molecule.residues]
    num_atoms_per_residue = [len(residue.atoms) for residue in molecule.residues]

    # Calculate SASA with ACCESS
    start_time = timeit.default_timer()
    res_atmarea, res_area, res_scarea = access.access1(num_residues, num_atoms, residue_ids, num_atoms_per_residue, atom_names, x, y, z, atom_radius, 1)
    elapsed = timeit.default_timer() - start_time

    return sum(res_atmarea), elapsed


if __name__ == "__main__":
    args = parse_command_line()

    with open(args.pdb_list) as input:
        for line in input:
            receptor, ligand = line.split()
            comp = receptor.replace('_rec.pdb', '.pdb')
            sasa_rec, time_rec = calculate_sasa(receptor)
            sasa_lig, time_lig = calculate_sasa(ligand)
            sasa_comp, time_comp = calculate_sasa(comp)
            print ','.join([str(sasa_rec), str(time_rec), str(sasa_lig), str(time_lig), str(sasa_comp), str(time_comp)])
