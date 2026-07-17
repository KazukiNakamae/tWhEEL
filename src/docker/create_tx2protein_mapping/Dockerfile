FROM python:3.12-slim

LABEL org.twheel.helper="create-tx2protein-mapping"

WORKDIR /app

COPY create_tx2protein_mapping.py /app/create_tx2protein_mapping.py

ENTRYPOINT ["python", "/app/create_tx2protein_mapping.py"]
