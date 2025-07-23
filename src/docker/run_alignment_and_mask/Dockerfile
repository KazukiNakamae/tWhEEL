FROM python:3.9-slim

# 必要なシステムパッケージをインストール
RUN apt-get update && apt-get install -y \
    mafft \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# 作業ディレクトリの作成
WORKDIR /app

# Pythonパッケージのインストール
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# アプリケーションファイルをコピー
COPY run_alignment_and_mask.py .

# 実行権限を付与
RUN chmod +x run_alignment_and_mask.py

# デフォルトコマンド
CMD ["python", "/app/run_alignment_and_mask.py"]
