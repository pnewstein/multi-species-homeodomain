BPC_IMGS = imgs/bpc_dots.svg imgs/bpc_umap.svg

CONDA_CMD = micromamba run -n mshdtf-env

all: $(BPC_IMGS)

$(BPC_IMGS): src/make_dotplots.py data/bpc.h5ad data/mouse_hdtf.csv
	$(CONDA_CMD) python src/make_dotplots.py bpc

data/mouse_hdtf.csv data/drosophila_hdtf.csv: src/get_hdtf.R
	mkdir -p data
	$(CONDA_CMD) R --no-save < src/get_hdtf.R

data/bpc.h5ad: src/make_bipolar_h5ad.py seq-data/clust_retinal_bipolar.txt seq-data/exp_matrix.txt
	$(CONDA_CMD) python src/make_bipolar_h5ad.py

seq-data/clust_retinal_bipolar.txt seq-data/exp_matrix.txt:
	cat seq-data/README.md >&2
