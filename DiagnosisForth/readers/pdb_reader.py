from typing import List, Union
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Mol
from tqdm import tqdm


def read_multiframe_pdb(file_path: Union[str, Path]) -> List[Mol]:
    """Read all frames from a multi-frame PDB file.
    
    Parameters
    ----------
    file_path : Union[str, Path]
        Path to the PDB file containing multiple frames/models
        
    Returns
    -------
    List[Mol]
        List of RDKit Mol objects, one for each frame in the PDB file
    """
    frame_contents: List[str] = []
    current_frame_lines: List[str] = []
    
    with open(file_path, 'r') as pdb_file:
        for line in tqdm(pdb_file):
            if line.startswith("MODEL"):
                current_frame_lines = []
            elif line.startswith("ENDMDL"):
                frame_contents.append(''.join(current_frame_lines))
            else:
                current_frame_lines.append(line)
                
    return [Chem.MolFromPDBBlock(frame_content) for frame_content in frame_contents]
