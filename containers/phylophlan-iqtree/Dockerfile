FROM mambaorg/micromamba:2.9.0

COPY --chown=$MAMBA_USER:$MAMBA_USER \
    modules/core/phylophlan/environment.yml /tmp/phylophlan-environment.yml

RUN micromamba install --yes --name base --file /tmp/phylophlan-environment.yml \
    && micromamba clean --all --yes

ENV PATH="/opt/conda/bin:${PATH}"
