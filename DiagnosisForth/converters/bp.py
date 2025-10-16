import sys
from typing import Union
from pathlib import Path
from rdkit import Chem
from Bio import PDB
from DiagnosisForth.utils.logger import log_sdf, log_pdb


def create_protein_ligand_complex(ligand_path: Union[str, Path],
                                  protein_path: Union[str, Path],
                                  output_path: Union[str, Path]) -> None:
    """Create a protein-ligand complex PDB file from separate protein and ligand files.
    
    Parameters
    ----------
    ligand_path : Union[str, Path]
        Path to the ligand file in SDF format
    protein_path : Union[str, Path]
        Path to the protein file in PDB format
    output_path : Union[str, Path]
        Path for the output complex PDB file
    """
    mol_supplier = Chem.SDMolSupplier(str(ligand_path), removeHs=False)
    ligand_mol = mol_supplier[0]
    if ligand_mol is None:
        print("Error: Couldn't read ligand.")
        sys.exit(1)
    
    try:
        smiles = Chem.MolToSmiles(ligand_mol)
        log_sdf(str(ligand_path), smiles, "create_protein_ligand_complex")
        log_pdb(str(protein_path), "create_protein_ligand_complex")
    except:
        pass

    pdb_parser = PDB.PDBParser(QUIET=True)
    protein_structure = pdb_parser.get_structure('protein', str(protein_path))

    ligand_pdb_block = Chem.MolToPDBBlock(ligand_mol)

    with open(output_path, 'w') as output_file:
        pdb_io = PDB.PDBIO()
        pdb_io.set_structure(protein_structure)
        pdb_io.save(output_file, write_end=False)
        output_file.write('TER\n')
        output_file.write(ligand_pdb_block)

    print("Done! Complex saved to", output_path)