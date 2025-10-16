import sys
from typing import Dict, List, Any, Tuple, Union, Optional

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

from plip.exchange.report import BindingSiteReport
from plip.structure.preparation import PDBComplex, LigandFinder
import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Mol
import tempfile




def analyze_protein_ligand_interactions_from_mol(pdb_mol: Mol) -> Dict[str, Dict[str, List[Any]]]:
    """Analyze protein-ligand interactions from an RDKit Mol object."""
    temp_file = tempfile.NamedTemporaryFile(
        mode='w+', delete=False, suffix='.pdb')
    Chem.MolToPDBFile(pdb_mol, temp_file.name)
    return analyze_protein_ligand_interactions(temp_file.name)




def analyze_protein_protein_interactions(pdb_file_path: Union[str, Path], chain1: str, chain2: str) -> Dict[str, Dict[str, List[Any]]]:
    """
    Analyze protein-protein interactions between specified chains using PLIP.

    Parameters
    ----------
    pdb_file_path : Union[str, Path]
        The PDB file of the complex.
    chain1 : str
        First protein chain identifier
    chain2 : str 
        Second protein chain identifier

    Returns
    -------
    dict :
        A dictionary of the binding sites and the interactions.
    """
    protein_complex = PDBComplex()
    protein_complex.load_pdb(str(pdb_file_path))
    ligand_finder = LigandFinder(protein_complex.protcomplex, protein_complex.altconf,
                      protein_complex.modres, protein_complex.covalent, protein_complex.Mapper)
    peptide_ligands = []
    for chain in [chain1, chain2]:
        peptide_ligands.append(ligand_finder.getpeptides(chain))
    for ligand in peptide_ligands:
        protein_complex.characterize_complex(ligand)

    interaction_sites = {}
    # loop over binding sites
    for key, site in sorted(protein_complex.interaction_sets.items()):
        # collect data about interactions
        binding_site = BindingSiteReport(site)
        # tuples of *_features and *_info will be converted to pandas DataFrame
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        # interactions is a dictionary which contains relevant information for each
        # of the possible interactions: hydrophobic, hbond, etc. in the considered
        # binding site. Each interaction contains a list with
        # 1. the features of that interaction, e.g. for hydrophobic:
        # ('RESNR', 'RESTYPE', ..., 'LIGCOO', 'PROTCOO')
        # 2. information for each of these features, e.g. for hydrophobic
        # (residue nb, residue type,..., ligand atom 3D coord., protein atom 3D coord.)
        interactions = {
            k: [getattr(binding_site, k + "_features")] +
            getattr(binding_site, k + "_info")
            for k in keys
        }
        interaction_sites[key] = interactions
    return interaction_sites


def analyze_protein_ligand_interactions(pdb_file_path: Union[str, Path], as_string: bool = False) -> Dict[str, Dict[str, List[Any]]]:
    """
    Analyze protein-ligand interactions using PLIP.

    Parameters
    ----------
    pdb_file_path : Union[str, Path]
        The PDB file of the complex.
    as_string : bool
        If True, treat pdb_file_path as PDB string content (not implemented)

    Returns
    -------
    dict :
        A dictionary of the binding sites and the interactions.
    """
    protein_complex = PDBComplex()
    protein_complex.load_pdb(str(pdb_file_path))  # load the pdb file
    for ligand in protein_complex.ligands:
        # find ligands and analyze interactions
        protein_complex.characterize_complex(ligand)
    interaction_sites = {}
    # loop over binding sites
    for key, site in sorted(protein_complex.interaction_sets.items()):
        # collect data about interactions
        binding_site = BindingSiteReport(site)
        # tuples of *_features and *_info will be converted to pandas DataFrame
        keys = (
            "hydrophobic",
            "hbond",
            "waterbridge",
            "saltbridge",
            "pistacking",
            "pication",
            "halogen",
            "metal",
        )
        # interactions is a dictionary which contains relevant information for each
        # of the possible interactions: hydrophobic, hbond, etc. in the considered
        # binding site. Each interaction contains a list with
        # 1. the features of that interaction, e.g. for hydrophobic:
        # ('RESNR', 'RESTYPE', ..., 'LIGCOO', 'PROTCOO')
        # 2. information for each of these features, e.g. for hydrophobic
        # (residue nb, residue type,..., ligand atom 3D coord., protein atom 3D coord.)
        interactions = {
            k: [getattr(binding_site, k + "_features")] +
            getattr(binding_site, k + "_info")
            for k in keys
        }
        interaction_sites[key] = interactions
    return interaction_sites


def protein_protein_interactions_to_dataframe(selected_site_interactions: Dict[str, List[Any]], 
                                             source_chain: str, 
                                             target_chain: str,
                                             interaction_type: str = "hbond") -> pd.DataFrame:
    """Convert protein-protein interaction data to pandas DataFrame."""
    interaction_df = interactions_to_dataframe(
        selected_site_interactions, interaction_type=interaction_type)
    interaction_df = interaction_df[(interaction_df['RESCHAIN'] == source_chain) & 
                                    (interaction_df['RESCHAIN_LIG'] == target_chain)]
    return interaction_df


def interactions_to_dataframe(binding_site: Dict[str, List[Any]], 
                             interaction_type: str = "hbond") -> pd.DataFrame:
    """
    Convert PLIP binding site data to pandas DataFrame.

    Parameters
    ----------
    binding_site : dict
        Precalculated interactions from PLIP for the selected site
    interaction_type : str
        The interaction type of interest (default set to hydrogen bond).

    Returns
    -------
    pd.DataFrame :
        DataFrame with information retrieved from PLIP.
    """

    # check if interaction type is valid:
    valid_types = [
        "hydrophobic",
        "hbond",
        "waterbridge",
        "saltbridge",
        "pistacking",
        "pication",
        "halogen",
        "metal",
    ]

    if interaction_type not in valid_types:
        print("!!! Wrong interaction type specified. Hbond is chosen by default!!!\n")
        interaction_type = "hbond"

    interaction_df = pd.DataFrame.from_records(
        # data is stored AFTER the column names
        binding_site[interaction_type][1:],
        # column names are always the first element
        columns=binding_site[interaction_type][0],
    )
    return interaction_df
