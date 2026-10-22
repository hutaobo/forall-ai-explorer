FROM python:3.11-slim
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY app.py core.py adjusted.py seekyoo_access.py sources.json prompt_completion.js ./
RUN mkdir -p /home/data && chown 1000:1000 /home/data
COPY public-data /data/forall
ENV FORALL_DATA_ROOT=/data/forall PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 NUMBA_NUM_THREADS=1 HOME=/tmp NUMBA_CACHE_DIR=/tmp/numba-cache XDG_CACHE_HOME=/tmp/.cache
ENV FORALL_SEEKYOO_AUTH=1 FORALL_PRIVATE_ROOT=/home/data/forall-private
EXPOSE 8501
USER 1000
ENTRYPOINT ["streamlit", "run", "app.py", "--server.address=0.0.0.0", "--server.port=8501", "--browser.gatherUsageStats=false"]
