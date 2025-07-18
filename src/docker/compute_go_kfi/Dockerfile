FROM python:3.9-slim

# 作業ディレクトリを設定
WORKDIR /app

# 必要なパッケージをインストール
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Pythonパッケージをインストール
RUN pip install pandas numpy scipy matplotlib seaborn argparse

# スクリプトをコピー
COPY compute_go_kfi.py /app/

# 実行権限を付与
RUN chmod +x /app/compute_go_kfi.py

# エントリーポイントを設定
ENTRYPOINT ["python", "/app/compute_go_kfi.py"]
