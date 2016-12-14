#!/usr/bin/env python

import os
import sys


if __name__ == "__main__":

    pdb_file = sys.argv[1]
    chains_rec = sys.argv[2]
    chains_lig = sys.argv[3]

    filename, file_extension = os.path.splitext(pdb_file)
    receptor_name = "%s_rec%s" % (filename, file_extension)
    ligand_name = "%s_lig%s" % (filename, file_extension)

    chains_rec = [chain.strip().upper() for chain in chains_rec.split(',')]
    chains_lig = [chain.strip().upper() for chain in chains_lig.split(',')]
    with open(pdb_file) as input:
        with open(receptor_name, 'w') as output_rec:
            with open(ligand_name, 'w') as output_lig:
                for line in input:
                    try:
                        if line:
                            if line[21] in chains_rec:
                                output_rec.write(line)
                            if line[21] in chains_lig:
                                output_lig.write(line)  
                    except:
                        pass
    print receptor_name, ligand_name
