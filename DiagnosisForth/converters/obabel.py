import subprocess
from typing import Union
from pathlib import Path
from DiagnosisForth.utils.logger import log_sdf, log_pdb


def convert_pdbqt_to_pdb(input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """Convert PDBQT format file to PDB format using OpenBabel."""
    command = f"obabel -ipdbqt {input_path} -opdb -O{output_path}"
    subprocess.run(command, shell=True)
    

def convert_pdbqt_to_sdf(input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """Convert PDBQT format file to SDF format using OpenBabel."""
    command = f"obabel -ipdbqt {input_path} -osdf -O{output_path}"
    subprocess.run(command, shell=True)
    
    try:
        from rdkit import Chem
        mol_supplier = Chem.SDMolSupplier(str(output_path), removeHs=False)
        if mol_supplier and len(mol_supplier) > 0:
            mol = mol_supplier[0]
            if mol:
                smiles = Chem.MolToSmiles(mol)
                log_sdf(str(output_path), smiles, "convert_pdbqt_to_sdf")
    except:
        pass