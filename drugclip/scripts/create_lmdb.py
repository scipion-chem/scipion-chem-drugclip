#!/usr/bin/env python3

import os
import argparse
import pickle
import lmdb
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser

def processMolecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[WARNING] Invalid SMILES: {smiles}")
        return None

    mol = Chem.AddHs(mol)
    status = AllChem.EmbedMolecule(mol, randomSeed=42)
    if status != 0:
        print(f"[WARNING] Failed to embed molecule: {smiles}")
        return None

    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())], dtype=np.float32)

    return {
        "atoms": atoms,
        "coordinates": coords,
        "mol": mol,
        "smi": smiles,
    }


def writeLmdb(entries, path):
    env = lmdb.open(
        path,
        map_size=1099511627776,
        subdir=False,
        lock=False
    )

    with env.begin(write=True) as txn:
        for i, entry in enumerate(tqdm(entries)):
            txn.put(str(i).encode(), pickle.dumps(entry))

    env.close()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--smiles-file", required=True,
                        help="Text file containing one SMILES per line")

    parser.add_argument("--pocket-files", required=True,
                        help="Comma-separated list of pocket PDB files")

    parser.add_argument("--output-dir", required=True,
                        help="Output directory for LMDB files")

    parser.add_argument("--max-pocket-atoms", type=int, default=256)

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Molecules
    molEntries = []
    with open(args.smiles_file) as f:
        for line in f:
            smiles = line.strip()
            entry = processMolecule(smiles)
            if entry:
                molEntries.append(entry)

    molLmdbPath = os.path.join(args.output_dir, "mols.lmdb")
    writeLmdb(molEntries, molLmdbPath)
    print(f"Molecules LMDB created: {molLmdbPath}")

    # Pockets
    pocketList = args.pocket_files.split(",")

    for pocketFile in pocketList:
        pocketFile = pocketFile.strip()
        pocketName = os.path.splitext(os.path.basename(pocketFile))[0]

        if pocketFile.endswith(".cif"):
            parserPocket = MMCIFParser(QUIET=True)
        else:  # pdb
            parserPocket = PDBParser(QUIET=True)

        structure = parserPocket.get_structure("pocket", pocketFile)

        pocketAtoms = []
        pocketCoords = []

        for atom in structure.get_atoms():
            pocketAtoms.append(atom.element)
            pocketCoords.append(list(atom.get_coord()))

        pocketAtoms = pocketAtoms[:args.max_pocket_atoms]
        pocketCoords = pocketCoords[:args.max_pocket_atoms]

        pocketEntry = {
            "pocket_atoms": pocketAtoms,
            "pocket_coordinates": np.array(pocketCoords, dtype=np.float32),
            "pocket": pocketName
        }
        pocketLmdbPath = os.path.join(
            args.output_dir, f"{pocketName}.lmdb"
        )

        writeLmdb([pocketEntry], pocketLmdbPath)  
        print(f"Pocket LMDB created: {pocketLmdbPath}")


if __name__ == "__main__":
    main()