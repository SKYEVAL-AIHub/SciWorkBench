FROM python:3.12-slim
ENV MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
WORKDIR /workspace
COPY requirements.txt /tmp/requirements.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt
COPY data/ /workspace/data/
RUN chmod -R a-w /workspace/data && mkdir -p /workspace/output
