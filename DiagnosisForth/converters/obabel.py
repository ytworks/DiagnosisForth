import sys
import subprocess
from typing import Union
from pathlib import Path

# Check consent before allowing module usage
from DiagnosisForth.utils.terms_consent import check_existing_consent, is_jupyter, request_consent_jupyter, request_consent

if not check_existing_consent():
    if is_jupyter():
        # Automatically show consent form in Jupyter
        consent_given = request_consent_jupyter()
        if not consent_given:
            raise ImportError("Terms of Service were declined. Cannot use this module.")
    else:
        # In terminal, request consent
        consent_given = request_consent()
        if not consent_given:
            sys.exit(1)

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