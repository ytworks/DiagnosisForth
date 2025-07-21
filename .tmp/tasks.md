# Task List

## Phase 1: 型ヒントの追加とカスタム型定義

### 1.1 カスタム型の定義
- [ ] kampo/types.pyファイルの作成
- [ ] InteractionPoint TypedDictの定義
- [ ] Interaction TypedDictの定義
- [ ] PharmacophoreFeature Literalの定義
- [ ] InteractionData型エイリアスの定義

### 1.2 converters モジュールの型ヒント追加
- [ ] converters/obabel.pyに型ヒントを追加
  - [ ] pdbqt2pdb関数
  - [ ] pdbqt2sdf関数
- [ ] converters/bp.pyに型ヒントを追加
  - [ ] getComplex関数

### 1.3 readers モジュールの型ヒント追加
- [ ] readers/pdb_reader.pyに型ヒントを追加
  - [ ] read_pdb_frames関数

### 1.4 interactions モジュールの型ヒント追加
- [ ] interactions/analyzer.pyに型ヒントを追加
  - [ ] retrieve_plip_interactions関数
  - [ ] retrieve_plip_interactions_for_ppi関数
  - [ ] create_df_from_binding_site関数

### 1.5 viewers モジュールの型ヒント追加
- [ ] viewers/utils.pyに型ヒントを追加
  - [ ] show_3d_interactions関数
  - [ ] show_3d_interactions_for_ppi関数
  - [ ] view_pharmacophore関数
  - [ ] extract_pharmacophore関数

## Phase 2: リネーミング

### 2.1 converters モジュールのリネーミング
- [ ] obabel.pyの関数リネーミング
  - [ ] pdbqt2pdb → convert_pdbqt_to_pdb
  - [ ] pdbqt2sdf → convert_pdbqt_to_sdf
  - [ ] 引数名の改善（filename → input_path, outname → output_path）
- [ ] bp.pyの関数リネーミング
  - [ ] getComplex → create_protein_ligand_complex
  - [ ] 引数名の改善（sdf → ligand_path, pdb → protein_path, out → output_path）

### 2.2 readers モジュールのリネーミング
- [ ] pdb_reader.pyの関数リネーミング
  - [ ] read_pdb_frames → read_multiframe_pdb
  - [ ] 引数名の改善（pdb_file → file_path）

### 2.3 interactions モジュールのリネーミング
- [ ] analyzer.pyの関数リネーミング
  - [ ] retrieve_plip_interactions → analyze_protein_ligand_interactions
  - [ ] retrieve_plip_interactions_for_ppi → analyze_protein_protein_interactions
  - [ ] create_df_from_binding_site → interactions_to_dataframe
  - [ ] 引数名の改善（pdb_file → pdb_file_path, c1/c2 → chain1/chain2）

### 2.4 viewers モジュールのリネーミング
- [ ] utils.pyの関数リネーミング
  - [ ] show_3d_interactions → visualize_protein_ligand_interactions
  - [ ] show_3d_interactions_for_ppi → visualize_protein_protein_interactions
  - [ ] view_pharmacophore → visualize_pharmacophore
  - [ ] extract_pharmacophore → extract_pharmacophore_features
  - [ ] 引数名の改善（highlight → highlighted_interactions）

### 2.5 __init__.pyファイルの更新
- [ ] 各モジュールの__init__.pyでエクスポートする関数名を更新
- [ ] kampo/__init__.pyでのインポートを更新

## Phase 3: テスト作成

### 3.1 テスト環境のセットアップ
- [ ] tests/ディレクトリの作成
- [ ] tests/__init__.pyの作成
- [ ] tests/conftest.pyの作成（共通フィクスチャ）
- [ ] pyproject.tomlにpytest関連の依存関係を追加
- [ ] pytest.iniまたはpyproject.tomlにpytest設定を追加

### 3.2 テストデータの準備
- [ ] tests/test_data/ディレクトリの作成
- [ ] サンプルPDB、SDF、PDBQTファイルの準備

### 3.3 converters モジュールのテスト
- [ ] tests/test_converters/ディレクトリの作成
- [ ] test_obabel.pyの作成
  - [ ] test_convert_pdbqt_to_pdb
  - [ ] test_convert_pdbqt_to_sdf
  - [ ] エラーケースのテスト
- [ ] test_bp.pyの作成
  - [ ] test_create_protein_ligand_complex
  - [ ] エラーケースのテスト

### 3.4 readers モジュールのテスト
- [ ] tests/test_readers/ディレクトリの作成
- [ ] test_pdb_reader.pyの作成
  - [ ] test_read_multiframe_pdb
  - [ ] 単一フレームPDBのテスト
  - [ ] 複数フレームPDBのテスト
  - [ ] エラーケースのテスト

### 3.5 interactions モジュールのテスト
- [ ] tests/test_interactions/ディレクトリの作成
- [ ] test_analyzer.pyの作成
  - [ ] test_analyze_protein_ligand_interactions
  - [ ] test_analyze_protein_protein_interactions
  - [ ] test_interactions_to_dataframe
  - [ ] 各種相互作用タイプのテスト

### 3.6 viewers モジュールのテスト
- [ ] tests/test_viewers/ディレクトリの作成
- [ ] test_utils.pyの作成
  - [ ] test_visualize_protein_ligand_interactions
  - [ ] test_visualize_protein_protein_interactions
  - [ ] test_visualize_pharmacophore
  - [ ] test_extract_pharmacophore_features

## Phase 4: サンプルコードの更新

### 4.1 Jupyter Notebookの更新
- [ ] examples/pco_example.ipynbの更新
  - [ ] 新しい関数名への置換
  - [ ] 型ヒントを活用したコード例
- [ ] examples/plip_example.ipynbの更新
  - [ ] 新しい関数名への置換
  - [ ] 型ヒントを活用したコード例
- [ ] examples/ppi.ipynbの更新
  - [ ] 新しい関数名への置換
  - [ ] 型ヒントを活用したコード例

### 4.2 動作確認
- [ ] 各サンプルノートブックの実行確認
- [ ] 出力結果の検証

## Phase 5: 品質保証

### 5.1 型チェック
- [ ] mypyの設定（pyproject.toml）
- [ ] 全モジュールでmypyを実行
- [ ] 型エラーの修正

### 5.2 コードフォーマット
- [ ] flake8またはruffでのリンティング
- [ ] blackでのフォーマット（必要に応じて）

### 5.3 テストカバレッジ
- [ ] pytest-covでカバレッジ測定
- [ ] カバレッジレポートの生成
- [ ] 重要な処理のカバレッジ確認

### 5.4 ドキュメント
- [ ] 各関数のdocstring確認・更新
- [ ] 型ヒントによる自己文書化の確認

## Phase 6: 最終確認

### 6.1 統合テスト
- [ ] エンドツーエンドのワークフロー確認
- [ ] 複数モジュールを組み合わせた処理の確認

### 6.2 後方互換性（オプション）
- [ ] 必要に応じて廃止予定のエイリアスを追加
- [ ] 移行ガイドの作成（必要に応じて）

## 推定作業時間

- Phase 1: 2-3時間
- Phase 2: 3-4時間
- Phase 3: 6-8時間
- Phase 4: 2-3時間
- Phase 5: 2-3時間
- Phase 6: 1-2時間

合計: 16-23時間