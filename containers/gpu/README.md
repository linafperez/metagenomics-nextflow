# Optional GPU images

GPU mode uses custom pinned images because the ordinary Biocontainers images
remain the safe CPU baseline. Build and publish the three contexts before use:

```bash
docker build -t metagenomics/comebin-gpu:1.0.4-cuda11.8 containers/comebin-gpu
docker build -t metagenomics/semibin-gpu:1.5.0-cuda11.8 containers/semibin-gpu
docker build -t metagenomics/vamb-gpu:5.0.4-cuda12.4 containers/vamb-gpu
```

The host must provide a compatible NVIDIA driver. Docker requires the NVIDIA
Container Toolkit; Apptainer/Singularity uses `--nv`. On SLURM the launcher also
requires the site-specific `--slurm-gpu-gres` value. Image builds and pulls are
not run by repository tests.
