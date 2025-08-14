# 要件定義: ライブラリ名変更（kampo → shishin）

## 1. 目的
- 現在のライブラリ名「kampo」を「shishin」に全面的に変更する

## 2. 変更範囲

### 2.1 ディレクトリ名
- `/kampo/` → `/shishin/`

### 2.2 プロジェクト設定
- `pyproject.toml`: パッケージ名の変更

### 2.3 ドキュメント・コメント
- 各Pythonファイルのdocstring内の「kampo」→「shishin」
- README.md内の全ての「kampo」参照

### 2.4 インポート文
- テストファイル内の全てのインポート文
- 例：`from kampo.xxx` → `from shishin.xxx`

### 2.5 Jupyterノートブック
- examples内の3つのノートブック：
  - pco_example.ipynb
  - ppi.ipynb
  - plip_example.ipynb
- インポート文、コメント、カーネル名の変更

### 2.6 その他
- .tmp内の既存ドキュメント（design.md、tasks.md）

## 3. 変更対象ファイル数
- Pythonファイル: 約10ファイル
- 設定ファイル: 1ファイル
- ノートブック: 3ファイル
- ドキュメント: 3ファイル

## 4. 制約事項
- 大文字小文字の区別を保持（Kampo → Shishin、kampo → shishin）
- 機能的な変更は行わない（名前の変更のみ）
- 全てのテストが変更後も正常に動作すること

## 5. 成功基準
- 全ての「kampo」参照が「shishin」に変更されていること
- プロジェクトが正常にインストール可能であること
- 全てのテストが合格すること
- ノートブックが正常に実行可能であること