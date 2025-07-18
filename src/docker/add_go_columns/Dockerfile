FROM python:3.9-slim

# 作業ディレクトリを設定
WORKDIR /app

# 必要なパッケージをインストール
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Pythonパッケージをインストール
RUN pip install pandas requests biopython

# スクリプトをコピー
COPY add_go_columns.py /app/

# 実行権限を付与
RUN chmod +x /app/add_go_columns.py

# エントリーポイントを設定
ENTRYPOINT ["python", "/app/add_go_columns.py"]
