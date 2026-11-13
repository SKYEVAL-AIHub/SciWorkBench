FROM python:3.12-slim
ENV MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
WORKDIR /app
COPY requirements.txt /tmp/requirements.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt
COPY data/ /app/data/
RUN chmod -R a-w /app/data && mkdir -p /app/output
