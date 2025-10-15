import subprocess
from typing import Union
from pathlib import Path


def convert_pdbqt_to_pdb(input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """Convert PDBQT format file to PDB format using OpenBabel."""
    command = f"obabel -ipdbqt {input_path} -opdb -O{output_path}"
    subprocess.run(command, shell=True)
    

def convert_pdbqt_to_sdf(input_path: Union[str, Path], output_path: Union[str, Path]) -> None:
    """Convert PDBQT format file to SDF format using OpenBabel."""
    command = f"obabel -ipdbqt {input_path} -osdf -O{output_path}"
    subprocess.run(command, shell=True)