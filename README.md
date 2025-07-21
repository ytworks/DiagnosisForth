# Kampo - Drug Discovery Toolkit

Kampo（漢方）は、創薬研究のためのPythonツールキットです。タンパク質-リガンド相互作用解析、ファーマコフォア特徴抽出、分子可視化などの機能を提供します。

## 主な機能

- **タンパク質-リガンド相互作用解析**: PLIPを使用した詳細な相互作用解析
- **タンパク質-タンパク質相互作用解析**: インターフェース残基の特定と相互作用の可視化
- **ファーマコフォア特徴抽出**: 分子構造からの薬理学的特徴の抽出
- **分子形式変換**: PDBQT形式からPDB/SDF形式への変換
- **3D分子可視化**: py3Dmolを使用したインタラクティブな3D表示

## インストール

### 前提条件

- Python 3.8以上
- OpenBabel（コマンドラインツール）

### インストール手順

1. **OpenBabelのインストール**
   ```bash
   # macOSの場合
   brew install open-babel
   
   # Ubuntuの場合
   sudo apt-get install openbabel
   
   # Condaを使用する場合
   conda install -c conda-forge openbabel
   ```

2. **Kampoのインストール**
   ```bash
   # リポジトリをクローン
   git clone https://github.com/yourusername/kampo.git
   cd kampo
   
   # pipでインストール
   pip install -e .
   ```

## クイックスタート

### タンパク質-リガンド相互作用解析

```python
from kampo.converters.obabel import convert_pdbqt_to_sdf
from kampo.converters.bp import create_protein_ligand_complex
from kampo.interactions.analyzer import analyze_protein_ligand_interactions
from kampo.viewers.utils import visualize_protein_ligand_interactions

# PDBQTファイルをSDF形式に変換
convert_pdbqt_to_sdf("ligand.pdbqt", "ligand.sdf")

# タンパク質-リガンド複合体を作成
create_protein_ligand_complex(
    ligand_path="ligand.sdf",
    protein_path="protein.pdb",
    output_path="complex.pdb"
)

# 相互作用を解析
interactions = analyze_protein_ligand_interactions("complex.pdb")

# 3Dで可視化
view, tables = visualize_protein_ligand_interactions("complex.pdb")
view.show()
```

### タンパク質-タンパク質相互作用解析

```python
from kampo.interactions.analyzer import analyze_protein_protein_interactions
from kampo.viewers.utils import visualize_protein_protein_interactions

# PPI解析を実行
interface_residues = analyze_protein_protein_interactions(
    "protein_complex.pdb",
    chain_a="A",
    chain_b="B"
)

# インターフェースを可視化
view = visualize_protein_protein_interactions(
    "protein_complex.pdb",
    chain_a="A", 
    chain_b="B",
    interface_residues=interface_residues
)
view.show()
```

### ファーマコフォア特徴抽出

```python
from kampo.interactions.analyzer import extract_pharmacophore_features
from kampo.viewers.utils import visualize_pharmacophores_from_mol

# 分子からファーマコフォア特徴を抽出
mol = Chem.MolFromSmiles("your_smiles_string")
features = extract_pharmacophore_features(mol)

# ファーマコフォアを可視化
view = visualize_pharmacophores_from_mol(mol, features)
view.show()
```

## サンプルノートブック

`examples/`ディレクトリに実践的な使用例があります：

- `plip_example.ipynb`: タンパク質-リガンド相互作用の詳細な解析
- `pco_example.ipynb`: ファーマコフォア特徴の抽出と可視化
- `ppi.ipynb`: タンパク質-タンパク質相互作用の解析

## 依存関係

主な依存関係：
- RDKit: 化学情報処理
- BioPython: タンパク質構造処理
- PLIP: タンパク質-リガンド相互作用解析
- py3Dmol: 3D分子可視化
- pandas: データ処理
- numpy: 数値計算
- matplotlib: グラフ作成

## 開発

### テストの実行

```bash
# すべてのテストを実行
pytest tests/

# カバレッジレポート付きで実行
pytest tests/ --cov=kampo
```

### コードスタイル

このプロジェクトはBlackとRuffを使用してコードフォーマットとリンティングを行っています：

```bash
# フォーマット
black kampo/

# リンティング
ruff check kampo/
```

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。

## 参考資料

- [TeachOpenCADD](https://projects.volkamerlab.org/teachopencadd/all_talktorials.html)
- [PLIP Documentation](https://github.com/pharmai/plip)
- [RDKit Documentation](https://www.rdkit.org/docs/)