all: imgs/bpc_dots.svg imgs/lamina_dots.svg imgs/AC.svg permute

CONDA_CMD = micromamba run -n mshdtf-env

BPC_OPT_INPUTS := $(wildcard seq-data/clust_retinal_bipolar.txt seq-data/exp_matrix.txt)
LAMINA_OPT_INPUTS := $(wildcard seq-data/Count)

permute: src/permute.py
	$(CONDA_CMD) python src/permute.py

data/adult_female_lamina_neurons.h5ad: src/make_lamina_h5ad.py $(LAMINA_OPT_INPUTS)
	@if [ -z "$(LAMINA_OPT_INPUTS)" ]; then \
		echo "Missing H5ad and .fastq.gz see README" >&2; \
		exit 1; \
	fi
	$(CONDA_CMD) python src/make_lamina_h5ad.py

imgs/bpc_dots.svg: src/make_dotplots.py data/panther_hdtf_families.txt data/bpc.h5ad
	mkdir -p imgs
	$(CONDA_CMD) python src/make_dotplots.py bpc

imgs/lamina_dots.svg: src/make_dotplots.py data/panther_hdtf_families.txt data/adult_female_lamina_neurons.h5ad
	mkdir -p imgs
	$(CONDA_CMD) python src/make_dotplots.py lamina

imgs/AC.svg imgs/HC.svg imgs/PR.svg imgs/RGC.svg: src/make_dotplots.py data/panther_hdtf_families.txt data/MRCA_full_normcounts.h5ad
	mkdir -p imgs
	$(CONDA_CMD) python src/make_dotplots.py retina

data/bpc.h5ad: src/make_bipolar_h5ad.py | $(BPC_OPT_INPUTS)
	@if [ -z "$(BPC_OPT_INPUTS)" ]; then \
		echo "Missing H5ad and matrix files see README" >&2; \
		exit 1; \
	fi
	mkdir -p data
	$(CONDA_CMD) python src/make_bipolar_h5ad.py

seq-data/clust_retinal_bipolar.txt seq-data/exp_matrix.txt: 
	cat seq-data/README.md >&2
	exit 1
