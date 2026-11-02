FROM mambaorg/micromamba:2.9.0
USER root
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ make git autoconf automake libtool pkg-config \
    && rm -rf /var/lib/apt/lists/*
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
WORKDIR /artifact
COPY --chown=$MAMBA_USER:$MAMBA_USER . /artifact
ENV PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
CMD ["python", "scripts/verify_submission.py"]
