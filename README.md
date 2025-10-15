# DiagnosisForth - Drug Discovery Toolkit

DiagnosisForth is a Python toolkit for drug discovery research, providing protein-ligand interaction analysis, pharmacophore feature extraction, and molecular visualization capabilities.

## Features

- **Protein-Ligand Interaction Analysis**: Detailed interaction analysis using PLIP
- **Protein-Protein Interaction Analysis**: Interface residue identification and interaction visualization
- **Pharmacophore Feature Extraction**: Extraction of pharmacological features from molecular structures
- **Molecular Format Conversion**: Conversion from PDBQT to PDB/SDF formats
- **3D Molecular Visualization**: Interactive 3D display using py3Dmol

## Installation

### Prerequisites

- Python 3.8 or higher
- OpenBabel (command-line tool)

### Installation Steps

1. **Install OpenBabel**
   ```bash
   # macOS
   brew install open-babel
   
   # Ubuntu
   sudo apt-get install openbabel
   
   # Using Conda
   conda install -c conda-forge openbabel
   ```

2. **Install DiagnosisForth**
   ```bash
   # Clone the repository
   git clone https://github.com/yourusername/DiagnosisForth.git
   cd DiagnosisForth
   
   # Install with pip
   pip install -e .
   ```

## Quick Start

### Protein-Ligand Interaction Analysis

```python
from DiagnosisForth.converters.obabel import convert_pdbqt_to_sdf
from DiagnosisForth.converters.bp import create_protein_ligand_complex
from DiagnosisForth.interactions.analyzer import analyze_protein_ligand_interactions
from DiagnosisForth.viewers.utils import visualize_protein_ligand_interactions

# Convert PDBQT file to SDF format
convert_pdbqt_to_sdf("ligand.pdbqt", "ligand.sdf")

# Create protein-ligand complex
create_protein_ligand_complex(
    ligand_path="ligand.sdf",
    protein_path="protein.pdb",
    output_path="complex.pdb"
)

# Analyze interactions
interactions = analyze_protein_ligand_interactions("complex.pdb")

# Visualize in 3D
view, tables = visualize_protein_ligand_interactions("complex.pdb")
view.show()
```

### Protein-Protein Interaction Analysis

```python
from DiagnosisForth.interactions.analyzer import analyze_protein_protein_interactions
from DiagnosisForth.viewers.utils import visualize_protein_protein_interactions

# Run PPI analysis
interface_residues = analyze_protein_protein_interactions(
    "protein_complex.pdb",
    chain_a="A",
    chain_b="B"
)

# Visualize interface
view = visualize_protein_protein_interactions(
    "protein_complex.pdb",
    chain_a="A", 
    chain_b="B",
    interface_residues=interface_residues
)
view.show()
```

### Pharmacophore Feature Extraction

```python
from DiagnosisForth.interactions.analyzer import extract_pharmacophore_features
from DiagnosisForth.viewers.utils import visualize_pharmacophores_from_mol

# Extract pharmacophore features from molecule
mol = Chem.MolFromSmiles("your_smiles_string")
features = extract_pharmacophore_features(mol)

# Visualize pharmacophores
view = visualize_pharmacophores_from_mol(mol, features)
view.show()
```

## Example Notebooks

Practical examples are available in the `examples/` directory:

- `plip_example.ipynb`: Detailed protein-ligand interaction analysis
- `pco_example.ipynb`: Pharmacophore feature extraction and visualization
- `ppi.ipynb`: Protein-protein interaction analysis

## Dependencies

Main dependencies:
- RDKit: Chemical informatics
- BioPython: Protein structure processing
- PLIP: Protein-ligand interaction profiling
- py3Dmol: 3D molecular visualization
- pandas: Data processing
- numpy: Numerical computing
- matplotlib: Plotting

## Development

### Running Tests

```bash
# Run all tests
pytest tests/

# Run with coverage report
pytest tests/ --cov=DiagnosisForth
```

### Code Style

This project uses Black and Ruff for code formatting and linting:

```bash
# Format code
black DiagnosisForth/

# Run linter
ruff check DiagnosisForth/
```

## License

This project is licensed under the GNU General Public License v2.0 - see the [LICENSE](LICENSE) file for details.

## References

- [TeachOpenCADD](https://projects.volkamerlab.org/teachopencadd/all_talktorials.html)
- [PLIP Documentation](https://github.com/pharmai/plip)
- [RDKit Documentation](https://www.rdkit.org/docs/)