# ワークフロー各ノードの説明

このディレクトリには、ワークフロー画像の各ノードに対応するPythonスクリプトが含まれています。

## ノード一覧

### Node 1: リガンドファイルダウンロード
- **ファイル**: `node_01_download_ligands.py`
- **入力**: なし（暗黙的）
- **出力**: `ligands.zip`
- **説明**: ZenodoからリガンドライブラリのZIPファイルをダウンロード

### Node 2: PDBファイルダウンロード
- **ファイル**: `node_02_download_pdb.py`
- **入力**: `pdb.id` (環境変数またはデフォルト値 "5Y7J")
- **出力**: `5Y7J.pdb`
- **説明**: RCSB PDBからタンパク質構造ファイルをダウンロード

### Node 3: リガンドファイル展開
- **ファイル**: `node_03_unpack_ligands.py`
- **入力**: `ligands.zip`
- **出力**: `constructed_library/` ディレクトリ（複数のSDFファイル）
- **説明**: ZIPファイルを展開してリガンドSDFファイルを取得

### Node 4: タンパク質入力確認
- **ファイル**: `node_04_protein_input.py`
- **入力**: `pdb.id` (環境変数またはデフォルト値)
- **出力**: `5Y7J.pdb` (確認のみ)
- **説明**: PDBファイルが存在することを確認

### Node 5: リガンド選択
- **ファイル**: `node_05_ligand_selection.py`
- **入力**: `constructed_library/` ディレクトリ内のSDFファイル群
- **出力**: `Ligands_select.sdf`
- **説明**: 選択されたリガンドを1つのSDFファイルにまとめる
- **注意**: デフォルトではテスト用に `clean_drug108*.sdf` パターンのファイルを選択

### Node 6: リガンドビュー
- **ファイル**: `node_06_ligand_view.py`
- **入力**: `Ligands_select.sdf`
- **出力**: `ligand.csv`
- **説明**: リガンドの情報（SMILES、分子量、LogPなど）をCSV形式で出力

### Node 7: Biopython - チェーン抽出
- **ファイル**: `node_07_extract_chains.py`
- **入力**: `5Y7J.pdb`
- **出力**: `5Y7J_chain.pdb`
- **説明**: PDBファイルからA鎖とB鎖、およびリファレンスリガンドを抽出

### Node 8: Biopython - リガンド中心識別
- **ファイル**: `node_08_ligand_center.py`
- **入力**: `Ligands_select.sdf`, `ligand.csv`, `5Y7J_chain.pdb`
- **出力**: `config.txt`
- **説明**: リガンドの中心座標を計算し、ドッキング設定ファイルを生成

### Node 9: OpenMM PDBFixer - タンパク質クリーンアップ
- **ファイル**: `node_09_clean_protein.py`
- **入力**: `5Y7J_chain.pdb`
- **出力**: `5Y7J_clean.pdb`, `5Y7J_AB_chains_fixed.pdb` (互換性のため)
- **説明**: PDBFixerを使用してタンパク質構造をクリーンアップ（欠損残基の追加、水素の追加など）

### Node 10: PDB2PQR - AMBER電荷適用
- **ファイル**: `node_10_apply_charges.py`
- **入力**: `5Y7J_clean.pdb`
- **出力**: `5Y7J_amber.pqr`, `5Y7J_amber.pdb`
- **説明**: PDB2PQRを使用してAMBER力場の電荷を適用

### Node 11: SMINA - インシリコスクリーニング
- **ファイル**: `node_11_smina_screening.py`
- **入力**: `config.txt`, `5Y7J_amber.pqr` (または `5Y7J_amber.pdb`), `Ligands_select.sdf`
- **出力**: `docking_results/` ディレクトリ（ドッキング結果ファイル群）
- **説明**: SMINAを使用して分子ドッキングを実行

### Node 12: レポート生成
- **ファイル**: `node_12_reporting.py`
- **入力**: `docking_results/` ディレクトリ内のドッキング結果
- **出力**: `results/docking_ranking.txt`, `results/*_docked.sdf` (トップ化合物)
- **説明**: ドッキング結果を解析してランキングを生成

## 実行順序

ワークフローは以下の順序で実行する必要があります：

```
1. Node 1 (リガンドダウンロード) ─┐
2. Node 2 (PDBダウンロード)      ─┤ 並列実行可能
                                   │
3. Node 3 (リガンド展開) ← Node 1
4. Node 4 (タンパク質入力確認) ← Node 2
5. Node 5 (リガンド選択) ← Node 3
6. Node 6 (リガンドビュー) ← Node 5
7. Node 7 (チェーン抽出) ← Node 2
                                   │
8. Node 8 (リガンド中心識別) ← Node 5, 6, 7
9. Node 9 (タンパク質クリーンアップ) ← Node 7
10. Node 10 (AMBER電荷適用) ← Node 9
                                  │
11. Node 11 (SMINAスクリーニング) ← Node 8, 10
12. Node 12 (レポート生成) ← Node 11
```

## 実行方法

### 個別実行

各ノードを個別に実行する場合：

```bash
python node_01_download_ligands.py
python node_02_download_pdb.py
python node_03_unpack_ligands.py
# ... 以下同様
```

### 一括実行

すべてのノードを順番に実行する場合：

```bash
# ステージ1: データ準備
python node_01_download_ligands.py
python node_02_download_pdb.py
python node_03_unpack_ligands.py
python node_04_protein_input.py

# ステージ2: リガンド処理
python node_05_ligand_selection.py
python node_06_ligand_view.py

# ステージ3: タンパク質処理
python node_07_extract_chains.py
python node_08_ligand_center.py
python node_09_clean_protein.py
python node_10_apply_charges.py

# ステージ4: ドッキング
python node_11_smina_screening.py

# ステージ5: レポート
python node_12_reporting.py
```

## 環境変数

- `PDB_ID`: PDB IDを指定（デフォルト: "5Y7J"）

```bash
export PDB_ID="5Y7J"
```

## 注意事項

1. **Node 5**: デフォルトではテスト用に `clean_drug108*.sdf` パターンのファイルを選択します。全ライブラリを処理する場合は、`node_05_ligand_selection.py` の `LIGAND_PATTERN` を変更してください。

2. **Node 11**: `Ligands_select.sdf` が存在しない場合、元のコードの動作（個別のSDFファイルを直接処理）にフォールバックします。

3. **依存関係**: 各ノードは前のノードの出力を入力として使用するため、順序を守って実行してください。

