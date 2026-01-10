# Workflow Functional Group Docker Environment

このディレクトリには、workflow-functional-groupの実行環境をDocker化するためのファイルが含まれています。

## ファイル構成

- `Dockerfile`: Ubuntu 24.04ベースのworkflow-functional-group環境のDockerイメージ定義
- `README.md`: このファイル

## 使用方法

### 1. Dockerイメージのビルド

```bash
# workflow-functional-groupディレクトリから実行
cd /home/shizuku/chiral/workflow-functional-group
docker build -t fegrow-env:latest -f 1-env/Dockerfile .
```

### 2. コンテナの実行

#### インタラクティブセッション
```bash
docker run -it --rm fegrow-env:latest
```

#### Jupyter Notebookサーバー
```bash
docker run -it --rm -p 8888:8888 fegrow-env:latest workflow-run jupyter notebook --ip=0.0.0.0 --port=8888 --no-browser --allow-root
```
ブラウザで `http://localhost:8888` にアクセス

#### 特定のコマンドを実行
```bash
docker run --rm fegrow-env:latest workflow-run python your_script.py
```

#### バックグラウンドで実行
```bash
docker run -d --name workflow-container fegrow-env:latest
```

#### 実行中のコンテナに接続
```bash
docker exec -it workflow-container /bin/bash
```

### 3. コンテナの管理

```bash
# 実行中のコンテナ一覧
docker ps

# すべてのコンテナ一覧（停止中も含む）
docker ps -a

# コンテナの停止
docker stop workflow-container

# コンテナの削除
docker rm workflow-container

# イメージの削除
docker rmi fegrow-env:latest
```

## 環境の詳細

このDocker環境には以下のものが含まれています：

- Ubuntu 24.04 ベース
- Miniconda + mamba
- workflow-functional-groupディレクトリ構造
- 必要なシステム依存関係

## 注意事項

- 初回ビルドには時間がかかります
- イメージサイズは約2.8GBです
- GPUを使用する場合は、`--gpus all` オプションを追加してください