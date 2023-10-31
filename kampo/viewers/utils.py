from ..interactions import analyzer
import pandas as pd
import py3Dmol


def show_3d_interactions(pdb_file, ligand_name='UNL'):
    """Show 3D interactions of a complex.

    Parameters
    ----------
    pdb_file : str
        Path to a PDB file of a complex.
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

    interactions_by_site = analyzer.retrieve_plip_interactions(pdb_file)
    index_of_selected_site = 0
    selected_site = list(interactions_by_site.keys())[index_of_selected_site]
    dfs = []
    complete_dfs = {}
    for bond_type in bond_types:
        df = analyzer.create_df_from_binding_site(
            interactions_by_site[selected_site],
            interaction_type=bond_type)
        complete_dfs[bond_type] = df
        df = df[["RESNR", "RESTYPE", "LIGCOO", "PROTCOO"]]
        df['BONDTYPE'] = bond_type
        dfs.append(df)
    df = pd.concat(dfs, axis=0)

    view = py3Dmol.view()
    view.addModel(open(pdb_file, 'r').read(), 'pdb')
    view.setStyle({"cartoon": {"color": "grey"}})
    LIG = [ligand_name]
    view.addStyle({'and': [{'resn': LIG}]},
                  {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
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
                          })
    view.setViewStyle({'style': 'outline', 'color': 'black', 'width': 0.1})

    view.zoomTo()
    return view, complete_dfs
