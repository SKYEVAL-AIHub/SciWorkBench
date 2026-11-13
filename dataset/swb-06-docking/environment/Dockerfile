FROM python:3.11-slim
ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg DEBIAN_FRONTEND=noninteractive
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends autodock-vina openbabel && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==2.2.6 pandas==2.3.3 scipy==1.16.3 matplotlib==3.10.8 rdkit==2025.3.3
COPY data /app/data
RUN chmod -R a-w /app/data && mkdir -p /app/output
