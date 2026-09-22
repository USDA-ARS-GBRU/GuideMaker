FROM condaforge/miniforge3:26.7.2-0

# 1. Create non-root user and ensure their home directory is generated (-m instead of -M)
RUN groupadd -g 10001 guidemaker && \
    useradd -u 10001 -g guidemaker -m -s /sbin/nologin guidemaker

WORKDIR /app

# 2. Copy only environment definition first
COPY environment.yml /app/environment.yml

# 3. Create a dedicated env (NOT base)
RUN mamba env create -p /env -f /app/environment.yml && \
    mamba clean --all -f -y

# Ensure conda activates correctly
ENV PATH=/env/bin:$PATH

# Force Streamlit to use its newly permitted home directory for config storage
ENV STREAMLIT_CONFIG_DIR=/home/guidemaker/.streamlit

# 4. Copy source code
COPY . /app

# 5. Install your package via pip inside conda env
RUN pip install --no-cache-dir --root-user-action=ignore .

# 6. Permissions (Added ownership of the generated home directory)
RUN chown -R guidemaker:guidemaker /app /env /home/guidemaker

EXPOSE 8501

# 7. Healthcheck (Swapped to 'python' to ensure it uses the activated 3.12 environment runner)
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8501/_stcore/health', timeout=10)" || exit 1

USER guidemaker

# 8. Streamlit runner
CMD ["streamlit", "run", "guidemaker/data/app.py", "--server.port=8501", "--server.address=0.0.0.0"]
