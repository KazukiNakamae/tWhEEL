# ─────────────── base ───────────────
FROM python:3.11-slim

# ---- install scientific stack -----
RUN pip install --no-cache-dir pandas goatools "matplotlib<3.11" squarify

# ---- working dir ----
WORKDIR /app

# ---- copy scripts ----
COPY count_go_parent_terms.py .
COPY go_slim_mapping.py .

ENTRYPOINT ["python", "count_go_parent_terms.py"]
