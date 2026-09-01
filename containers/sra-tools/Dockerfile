FROM mambaorg/micromamba:2.3.3

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes

ENV PATH=/opt/conda/bin:${PATH}

RUN prefetch --version 2>&1 | grep -F '3.4.1' \
    && fasterq-dump --version 2>&1 | grep -F '3.4.1' \
    && pigz --version 2>&1 | grep -F '2.8'
