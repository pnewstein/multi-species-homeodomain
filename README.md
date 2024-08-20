## Get sequencing data
### BPC data
- H5ad files can be found at ...
- move both h5ad files to seq-data directory
### Generating H5AD files
- Download metadata and expression matrix from the [Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP3/retinal-bipolar-neuron-drop-seq#/)
- move both files into the seq-data directory

## INSTALATION
- Install [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
- Clone the reposistory
```bash
git clone <>
```
- Create a new environment with installed dependencies
```bash
micromamba env create -f environment.yml
```

### Troubleshooting
- May fail to install on osx-arm64 
    - see [issue](https://github.com/conda-forge/r-xml-feedstock/issues/47)
    - instead use following commands:
    ```bash
micromamba env create -f backup-environment.yml
micromamba run -n mshdtf-env R --no-save < backup-biomart-installer.R
    ```
