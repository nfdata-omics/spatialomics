FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="cytassist-microscopy-registration"
LABEL org.opencontainers.image.description="Python image-registration environment for CytAssist and microscopy image alignment"
LABEL org.opencontainers.image.source="https://github.com/nfdata-omics/spatialomics"

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes

ENV MPLBACKEND=Agg

CMD ["python", "--version"]
