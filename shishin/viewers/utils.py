from typing import Dict, List, Any, Tuple, Union, Optional
from ..interactions import analyzer
import pandas as pd
import py3Dmol
import tempfile
from rdkit import Chem
from rdkit.Chem import Mol
from rdkit.Chem import ChemicalFeatures
from rdkit.Chem import AllChem
from rdkit import RDConfig
import os
from pathlib import Path


pharmacophore_colors: Dict[str, Tuple[int, int, int]] = {
    "Donor": (0, 255, 0),        # 緑色
    "Acceptor": (255, 0, 0),     # 赤色
    "NegIonizable": (0, 0, 255),  # 青色
    "PosIonizable": (255, 255, 0),  # 黄色
    "Aromatic": (128, 0, 128),   # 紫色
    "ZnBinder": (0, 255, 255),   # シアン色
    "Hydrophobe": (255, 165, 0),  # オレンジ色
    "LumpedHydrophobe": (255, 105, 180)  # ピンク色
}


def visualize_protein_ligand_interactions_from_mol(pdb_mol: Mol, ligand_name: str = 'UNL') -> Tuple[py3Dmol.view, Dict[str, pd.DataFrame]]:
    """Visualize protein-ligand interactions from an RDKit Mol object."""
    temp_file = tempfile.NamedTemporaryFile(
        mode='w+', delete=False, suffix='.pdb')
    Chem.MolToPDBFile(pdb_mol, temp_file.name)
    return visualize_protein_ligand_interactions(temp_file.name)


def visualize_pharmacophore_from_mol(pdb_mol: Mol,
                                    smiles: str = 'O=C(O)CCCN1CC(Oc2c1cccc2NC(=O)c1ccc(cc1)OCCCCCc1ccccc1)C(=O)O',
                                    ligand_name: str = "LIG") -> Tuple[py3Dmol.view, Dict[str, pd.DataFrame], pd.DataFrame]:
    """Visualize pharmacophore features from an RDKit Mol object."""
    view, complete_dfs = visualize_protein_ligand_interactions_from_mol(pdb_mol, ligand_name)
    temp_file = tempfile.NamedTemporaryFile(
        mode='w+', delete=False, suffix='.pdb')
    Chem.MolToPDBFile(pdb_mol, temp_file.name)
    ligand_mol = extract_ligand_from_pdb(temp_file.name,
                                         smiles=smiles,
                                         residue_name=ligand_name)
    pharmacophore_features = extract_pharmacophore_features(ligand_mol=ligand_mol)
    view = add_pharmacophore_spheres(view, pharmacophore_features)
    features_df = pd.DataFrame(
        pharmacophore_features, columns=['Type', 'Detail', 'AtomNumber', 'X', 'Y', 'Z'])

    return view, complete_dfs, features_df


def add_pharmacophore_spheres(view: py3Dmol.view, pharmacophore_features: List[List[Any]]) -> py3Dmol.view:
    """Add pharmacophore feature spheres to the 3D view."""
    for feature in pharmacophore_features:
        color = pharmacophore_colors[feature[0]]
        view.addSphere({'center': {'x': feature[3], 'y': feature[4], 'z': feature[5]},
                        'radius': 0.5,
                        'color': f'rgb({color[0]},{color[1]},{color[2]})',
                        "hoverable": True,
                        "hover_callback": '''function(atom,viewer,event,container) {
                                                if(!this.label) {
                                                    this.label = viewer.addLabel("%s = (%s, %s, %s)",{position: this, backgroundColor: 'mintcream', fontColor:'black'});
                                            }}''' % (feature[0], feature[3], feature[4], feature[5]),
                        "unhover_callback": '''function(atom,viewer) { 
                                    if(this.label) {
                                        viewer.removeLabel(this.label);
                                        delete this.label;
                                    }
                                    }'''})
    return view


def visualize_pharmacophore(pdb_file: Union[str, Path],
                           smiles: str = 'O=C(O)CCCN1CC(Oc2c1cccc2NC(=O)c1ccc(cc1)OCCCCCc1ccccc1)C(=O)O',
                           ligand_name: str = "LIG") -> Tuple[py3Dmol.view, Dict[str, pd.DataFrame], pd.DataFrame]:
    """Visualize pharmacophore features of a ligand."""
    view, complete_dfs = visualize_protein_ligand_interactions(pdb_file, ligand_name)
    ligand_mol = extract_ligand_from_pdb(pdb_file,
                                         smiles=smiles,
                                         residue_name=ligand_name)
    pharmacophore_features = extract_pharmacophore_features(ligand_mol=ligand_mol)
    view = add_pharmacophore_spheres(view, pharmacophore_features)
    features_df = pd.DataFrame(
        pharmacophore_features, columns=['Type', 'Detail', 'AtomNumber', 'X', 'Y', 'Z'])
    return view, complete_dfs, features_df


def extract_ligand_from_pdb(pdb_file: Union[str, Path], 
                            smiles: str = 'O=C(O)CCCN1CC(Oc2c1cccc2NC(=O)c1ccc(cc1)OCCCCCc1ccccc1)C(=O)O',
                            residue_name: str = "LIG") -> Mol:
    mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
    ligand_atoms = [atom.GetIdx() for atom in mol.GetAtoms(
    ) if atom.GetPDBResidueInfo().GetResidueName() == residue_name]

    # リガンド原子のみを含むサブモレキュールを作成
    ligand = Chem.RWMol(mol)
    ligand_atoms_set = set(ligand_atoms)
    for atom_idx in reversed(range(ligand.GetNumAtoms())):
        if atom_idx not in ligand_atoms_set:
            ligand.RemoveAtom(atom_idx)
    lig = ligand.GetMol()
    reference_mol = Chem.MolFromSmiles(smiles)
    ligand_mol = AllChem.AssignBondOrdersFromTemplate(
        reference_mol, lig)
    ligand_mol.AddConformer(lig.GetConformer(0))
    return ligand_mol


def visualize_protein_protein_interactions(pdb_file_path: Union[str, Path], 
                                          chain1: str, 
                                          chain2: str, 
                                          highlighted_interactions: Optional[List[str]] = None,
                                          width: int = 900,
                                          height: int = 600) -> Tuple[py3Dmol.view, Dict[str, pd.DataFrame]]:
    bond_types = ["hydrophobic",
                  "hbond",
                  "waterbridge",
                  "saltbridge",
                  "pistacking",
                  "pication",
                  "halogen"]
    colors = {"hydrophobic": "red",
              "hbond": "blue",
              "waterbridge": "blueCarbon",
              "saltbridge": "white",
              "pistacking": "yellow",
              "pication": "orange",
              "halogen": "magenta"}
    interactions_by_site = analyzer.analyze_protein_protein_interactions(
        pdb_file_path, chain1, chain2)
    index_of_selected_site = 0
    selected_site = list(interactions_by_site.keys())[index_of_selected_site]
    dfs = []
    complete_dfs = {}
    for bond_type in bond_types:
        df = analyzer.protein_protein_interactions_to_dataframe(
            interactions_by_site[selected_site],
            source_chain=chain1,
            target_chain=chain2,
            interaction_type=bond_type)
        complete_dfs[bond_type] = df
        df = df[["RESNR", "RESTYPE", "LIGCOO",
                 "PROTCOO", 'RESNR_LIG', "RESTYPE_LIG"]]
        df['BONDTYPE'] = bond_type
        dfs.append(df)
    df = pd.concat(dfs, axis=0)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(open(pdb_file_path, 'r').read(), 'pdb')
    view.setStyle({"cartoon": {"color": "grey"}})
    view = add_hover_callbacks(view)
    for idx, rows in df.iterrows():
        view.setStyle({'resi': rows["RESNR"]},
                      {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.2}})
        view.setStyle({'resi': rows["RESNR_LIG"]},
                      {'stick': {'colorscheme': 'whiteCarbon', 'radius': 0.2}})
        sx, sy, sz = rows["LIGCOO"]
        ex, ey, ez = rows["PROTCOO"]
        view.addCylinder({"start": dict(x=sx, y=sy, z=sz),
                          "end":   dict(x=ex, y=ey, z=ez),
                          "color": colors[rows["BONDTYPE"]],
                          "radius": .15,
                          "dashed": True,
                          "fromCap": 1,
                          "toCap": 1,
                          "hoverable": True,
                          "hover_callback": '''function(atom,viewer,event,container) {
                                                if(!this.label) {
                                                    this.label = viewer.addLabel("%s",{position: this, backgroundColor: 'mintcream', fontColor:'black'});
                                            }}''' % rows["BONDTYPE"],
                          "unhover_callback": '''function(atom,viewer) { 
                                    if(this.label) {
                                        viewer.removeLabel(this.label);
                                        delete this.label;
                                    }
                                    }'''
                          })
    view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.1})

    view.zoomTo()
    return view, complete_dfs


def visualize_protein_ligand_interactions(pdb_file_path: Union[str, Path], 
                                         interaction_data: Optional[Dict[str, List[Any]]] = None,
                                         highlighted_interactions: Optional[List[str]] = None,
                                         width: int = 900,
                                         height: int = 600) -> Tuple[py3Dmol.view, Dict[str, pd.DataFrame]]:
    """Visualize protein-ligand interactions in 3D.

    Parameters
    ----------
    pdb_file_path : Union[str, Path]
        Path to a PDB file of a complex.
    interaction_data : Optional[Dict[str, List[Any]]] 
        Pre-calculated interaction data (if None, will be calculated)
    highlighted_interactions : Optional[List[str]]
        List of interaction types to highlight
    width : int
        Width of the viewer
    height : int
        Height of the viewer
    """
    bond_types = ["hydrophobic",
                  "hbond",
                  "waterbridge",
                  "saltbridge",
                  "pistacking",
                  "pication",
                  "halogen"]
    colors = {"hydrophobic": "red",
              "hbond": "blue",
              "waterbridge": "blueCarbon",
              "saltbridge": "white",
              "pistacking": "yellow",
              "pication": "orange",
              "halogen": "magenta"}

    if interaction_data is None:
        interactions_by_site = analyzer.analyze_protein_ligand_interactions(pdb_file_path)
    else:
        interactions_by_site = interaction_data
    index_of_selected_site = 0
    selected_site = list(interactions_by_site.keys())[index_of_selected_site]
    dfs = []
    complete_dfs = {}
    for bond_type in bond_types:
        df = analyzer.interactions_to_dataframe(
            interactions_by_site[selected_site],
            interaction_type=bond_type)
        complete_dfs[bond_type] = df
        df = df[["RESNR", "RESTYPE", "LIGCOO", "PROTCOO"]]
        df['BONDTYPE'] = bond_type
        dfs.append(df)
    df = pd.concat(dfs, axis=0)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(open(pdb_file_path, 'r').read(), 'pdb')
    view.setStyle({"cartoon": {"color": "grey"}})
    view = add_hover_callbacks(view)
    # Note: ligand_name detection would need to be implemented
    ligand_residues = ['UNL', 'LIG']  # Common ligand names
    view.addStyle({'and': [{'resn': ligand_residues}]},
                  {'stick': {'colorscheme': 'magentaCarbon', 'radius': 0.3}})
    view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.1})
    for idx, rows in df.iterrows():
        view.setStyle({'resi': rows["RESNR"]},
                      {'stick': {'colorscheme': 'whiteCarbon', 'radius': 0.2}})
        sx, sy, sz = rows["LIGCOO"]
        ex, ey, ez = rows["PROTCOO"]
        view.addCylinder({"start": dict(x=sx, y=sy, z=sz),
                          "end":   dict(x=ex, y=ey, z=ez),
                          "color": colors[rows["BONDTYPE"]],
                          "radius": .15,
                          "dashed": True,
                          "fromCap": 1,
                          "toCap": 1,
                          "hoverable": True,
                          "hover_callback": '''function(atom,viewer,event,container) {
                                                if(!this.label) {
                                                    this.label = viewer.addLabel("%s",{position: this, backgroundColor: 'mintcream', fontColor:'black'});
                                            }}''' % rows["BONDTYPE"],
                          "unhover_callback": '''function(atom,viewer) { 
                                    if(this.label) {
                                        viewer.removeLabel(this.label);
                                        delete this.label;
                                    }
                                    }'''
                          })
    view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.1})

    view.zoomTo()
    return view, complete_dfs


def extract_pharmacophore_features(ligand_mol: Mol) -> List[List[Any]]:
    """Extract pharmacophore features from a ligand molecule."""
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    pharmacophore_features = []

    features = factory.GetFeaturesForMol(ligand_mol)
    for feature in features:
        pharmacophore_features.append([
            feature.GetFamily(), 
            feature.GetType(), 
            feature.GetAtomIds(),
            feature.GetPos()[0], 
            feature.GetPos()[1], 
            feature.GetPos()[2]
        ])
    return pharmacophore_features


def add_hover_callbacks(view: py3Dmol.view) -> py3Dmol.view:
    """Add hover callbacks to atoms in the 3D view."""
    view.setHoverable({}, True, '''function(atom,viewer,event,container) {
                   if(!atom.label) {
                    atom.label = viewer.addLabel(atom.resn+":"+atom.resi,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                   }}''',
                      '''function(atom,viewer) { 
                   if(atom.label) {
                    viewer.removeLabel(atom.label);
                    delete atom.label;
                   }
                }''')

    return view
