FROM mambaorg/micromamba:2.3.3
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes
ENV PATH=/opt/conda/bin:${PATH}
RUN python3 -c "import torch; assert torch.__version__.split('+', 1)[0] == '2.1.2'; assert torch.version.cuda == '11.8'"
RUN run_comebin.sh --help >/dev/null
