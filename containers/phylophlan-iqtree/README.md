# PhyloPhlAn and IQ-TREE container

This build definition creates the exact combined runtime required by the
phylogenomics module: PhyloPhlAn 3.1.1 and IQ-TREE 3.0.1.

Build from the repository root and publish the image to a registry accessible
from every execution node:

```bash
docker build \
    --file containers/phylophlan-iqtree/Dockerfile \
    --tag metagenomics/phylophlan-iqtree:3.1.1-3.0.1 \
    .

```

That tag is the pipeline default for local Docker. For Apptainer, Singularity,
or distributed Docker execution, retag and publish it to an accessible OCI
registry, then pass that image to the pipeline:

```bash
--phylophlan_container registry.example.org/phylophlan-iqtree:3.1.1-3.0.1
```

Apptainer and Singularity can consume the same OCI image through Nextflow.
