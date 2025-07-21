# Design Document

## System Architecture Overview

kampoは創薬研究向けのPythonライブラリで、以下のモジュール構成となっています：

```
kampo/
├── converters/    # ファイル形式変換
├── interactions/  # 相互作用解析
├── readers/       # ファイル読み込み
└── viewers/       # 3D可視化
```

## Component Design

### 1. converters モジュール

#### 1.1 obabel.py
**現在の問題点:**
- 関数名が簡潔すぎる（`pdbqt2pdb`, `pdbqt2sdf`）
- 型ヒントなし
- エラーハンドリングが不十分

**リファクタリング設計:**
```python
from typing import Optional, Union
from pathlib import Path

def convert_pdbqt_to_pdb(
    input_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None
) -> Path:
    """Convert PDBQT format file to PDB format."""
    pass

def convert_pdbqt_to_sdf(
    input_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None
) -> Path:
    """Convert PDBQT format file to SDF format."""
    pass
```

#### 1.2 bp.py
**現在の問題点:**
- 関数名が不明瞭（`getComplex`）
- camelCaseを使用（PEP 8違反）

**リファクタリング設計:**
```python
def create_protein_ligand_complex(
    ligand_path: Union[str, Path],
    protein_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None
) -> Path:
    """Create a protein-ligand complex PDB file."""
    pass
```

### 2. readers モジュール

#### 2.1 pdb_reader.py
**現在の問題点:**
- 関数名が汎用的すぎる（`read_pdb_frames`）

**リファクタリング設計:**
```python
from typing import List
from rdkit.Chem import Mol

def read_multiframe_pdb(
    file_path: Union[str, Path]
) -> List[Mol]:
    """Read all frames from a multi-frame PDB file."""
    pass
```

### 3. interactions モジュール

#### 3.1 analyzer.py
**現在の問題点:**
- 長い関数名だが一貫性がない
- 返り値の型が不明確
- DataFrameの構造が文書化されていない

**リファクタリング設計:**
```python
from typing import Dict, List, Any, Tuple
import pandas as pd
from plip.structure.preparation import PDBComplex

InteractionData = Dict[str, List[Dict[str, Any]]]

def analyze_protein_ligand_interactions(
    pdb_file_path: Union[str, Path],
    as_string: bool = False
) -> Tuple[PDBComplex, InteractionData]:
    """Analyze protein-ligand interactions using PLIP."""
    pass

def analyze_protein_protein_interactions(
    pdb_file_path: Union[str, Path],
    chain1: str,
    chain2: str
) -> Tuple[PDBComplex, InteractionData]:
    """Analyze protein-protein interactions between specified chains."""
    pass

def interactions_to_dataframe(
    binding_site: Any,
    interaction_type: str = "ligand"
) -> pd.DataFrame:
    """Convert PLIP binding site data to pandas DataFrame."""
    pass
```

### 4. viewers モジュール

#### 4.1 utils.py
**現在の問題点:**
- 関数名に冗長性（`show_3d_`プレフィックス）
- パラメータ名が不明瞭（`highlight`）

**リファクタリング設計:**
```python
from typing import Optional, List, Dict, Any
import py3Dmol

def visualize_protein_ligand_interactions(
    pdb_file_path: Union[str, Path],
    interaction_data: InteractionData,
    highlighted_interactions: Optional[List[str]] = None,
    width: int = 900,
    height: int = 600
) -> py3Dmol.view:
    """Visualize protein-ligand interactions in 3D."""
    pass

def visualize_protein_protein_interactions(
    pdb_file_path: Union[str, Path],
    interaction_data: InteractionData,
    chain1: str,
    chain2: str,
    highlighted_interactions: Optional[List[str]] = None,
    width: int = 900,
    height: int = 600
) -> py3Dmol.view:
    """Visualize protein-protein interactions in 3D."""
    pass

def visualize_pharmacophore(
    ligand_mol: Union[str, Path, Mol],
    features: Optional[List[str]] = None,
    width: int = 400,
    height: int = 400
) -> py3Dmol.view:
    """Visualize pharmacophore features of a ligand."""
    pass

def extract_pharmacophore_features(
    ligand_mol: Union[str, Path, Mol]
) -> Dict[str, List[Tuple[float, float, float]]]:
    """Extract pharmacophore features from a ligand molecule."""
    pass
```

## Data Structures

### Custom Types
```python
# types.py
from typing import TypedDict, List, Tuple, Literal

class InteractionPoint(TypedDict):
    residue: str
    chain: str
    position: int
    atom: str
    coordinates: Tuple[float, float, float]

class Interaction(TypedDict):
    type: Literal[
        "hydrophobic", "hydrogen_bond", "water_bridge",
        "salt_bridge", "pi_stacking", "pi_cation", "halogen"
    ]
    protein: InteractionPoint
    ligand: InteractionPoint
    distance: float
    angle: Optional[float]

PharmacophoreFeature = Literal[
    "Donor", "Acceptor", "Hydrophobe", 
    "PosIonizable", "NegIonizable", "Aromatic"
]
```

## Testing Strategy

### Test Structure
```
tests/
├── __init__.py
├── conftest.py              # pytest fixtures
├── test_converters/
│   ├── __init__.py
│   ├── test_obabel.py
│   └── test_bp.py
├── test_readers/
│   ├── __init__.py
│   └── test_pdb_reader.py
├── test_interactions/
│   ├── __init__.py
│   └── test_analyzer.py
└── test_viewers/
    ├── __init__.py
    └── test_utils.py
```

### Test Categories
1. **Unit Tests**: 各関数の個別テスト
2. **Integration Tests**: モジュール間の連携テスト
3. **Edge Cases**: 異常系・境界値テスト

## Migration Plan

### Phase 1: 型ヒントの追加
1. カスタム型の定義（types.py）
2. 各モジュールへの型ヒント追加

### Phase 2: リネーミング
1. 関数名の変更（snake_case、意味のある名前）
2. 変数名の変更
3. 引数名の変更

### Phase 3: テスト作成
1. テスト環境のセットアップ
2. 各モジュールのユニットテスト作成
3. 統合テストの作成

### Phase 4: サンプルコードの更新
1. examples内のJupyterノートブックの更新
2. 新しい関数名への置換
3. 動作確認

## Dependencies

### Required Updates
- pyproject.tomlに依存関係を明記
- 型チェック用にmypyを追加
- テスト用にpytest関連パッケージを追加

### Development Dependencies
```toml
[tool.poetry.dev-dependencies]
pytest = "^7.0"
pytest-cov = "^4.0"
mypy = "^1.0"
types-setuptools = "*"
```

## Naming Convention Guidelines

### Functions
- 動詞で始まる（analyze_, convert_, create_, visualize_）
- 処理内容を明確に表現
- snake_case使用

### Variables
- 意味のある名前（`df` → `interaction_dataframe`）
- 略語を避ける（`pdb` → `pdb_file_path`）

### Constants
- 大文字とアンダースコア（`DEFAULT_WIDTH`）

## Backwards Compatibility

必要に応じて、以下の方法で後方互換性を提供：
```python
# Deprecated aliases
pdbqt2pdb = convert_pdbqt_to_pdb  # 廃止予定の警告付き
```