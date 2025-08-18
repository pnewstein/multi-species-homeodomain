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

# Source of panther_hdtf_family.txt
download from
https://www.pantherdb.org/panther/familyList.do?searchType=basic&fieldName=all&listType=6&fieldValue=PC00119

ensure that all items are displayed
click send list to file
save
take first word from each line that does not have a subfamily
awk '{print $1}' pantherGeneList.txt | grep -v : > panther_hdtf_families.txt

