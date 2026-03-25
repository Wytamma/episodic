# Running Episodic on a SLURM cluster

> Episodic + SLURM + GPU go brrr. 

Episodic has built-in support for SLURM clusters via Snakemake profiles, allowing you to run this BEAST analyses on high-performance computing resources. This guide shows how to configure and execute the workflow on a SLURM cluster, including GPU-accelerated BEAST runs.

## SLURM profiles

Episodic includes two Snakemake profiles for SLURM:
- `slurm`: for CPU-based BEAST runs.
- `slurm-gpu`: for GPU-accelerated BEAST runs using BEAGLE

To use these profiles, specify `--profile slurm` or `--profile slurm-gpu` when running the workflow. The profiles will automatically submit BEAST jobs to the SLURM scheduler according to the resource requirements defined in the workflow configuration.

You can use the episodic profile command to generate a profile template and customize it for your cluster:

```console
episodic profile show slurm-gpu > my-slurm-gpu/config.yaml
```

Then run `episodic run --profile my-slurm-gpu` to execute the workflow with your custom SLURM settings.

!!! info

    You will likely need to customise the gpu partition name as it defaults to `gpu-a100` in the profile. You can also adjust other resource settings such as memory and runtime limits as needed (see below).

## Overriding SLURM settings via the cli with snakemake flags

You can use `--set-resources` to override SLURM resource settings on the command line. For example, to request 16 GB of memory and a runtime of 4 hours for BEAST runs:

```console
episodic run --profile slurm --set-resources beast:mem_mb=16G --set-resources beast:runtime=04:00:00
```

!!! note
    
    The marginal likelihood estimation rule (`mle`) also has resource settings that can be overridden in the same way, for example:
    
    ```console
    episodic run --profile slurm --set-resources mle:mem_mb=16G --set-resources mle:runtime=04:00:00 --set-resources run_beast:partition=gpu-a100
    ```

## Optimising GPU runs

When running beast on GPU nodes its often beneficial to run multiple BEAST runs per node to maximise GPU utilisation. You can achieve this by using the [groups feature of Snakemake](https://snakemake.readthedocs.io/en/stable/executing/grouping.html) to group multiple BEAST runs together and submit them as a single job. For example, to group BEAST runs in batches of 4:

```console
episodic run --profile slurm-gpu --groups beast=gpuGroup --group-components gpuGroup=4
```