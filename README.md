![Echinus esculentus](img/cover.png)

# echinus

Get a gene expression table from counts for equivalence classes (ECs -> GE table)

Echinus provides a fast & easy way to summarise transcript compatibility counts to gene abundances. Counts for equivalence classes corresponding to one gene are summed up in order to create a cells x genes matrix commonly used for downstream analysis.

Please address documentation for [kallisto](https://pachterlab.github.io/kallisto/manual) and [bustools](https://bustools.github.io/manual) for details on getting the TCCs for your data. In particular, check out the possibility to use `kallisto bus` command to pseudo-align reads while handling your barcode-UMI design via `-x` argument.

# Usage

Echinus can be used on TCCs that are produced by `kallisto pseudo` command:

```
echinus count --tcc matrix.tcc.mtx --ecmap matrix.ec --txnames transcripts.txt --genemap transcripts_to_genes.txt --output echinus_dir 
```

or on the output of the `kallisto|bustools` pipeline after correcting the barcodes and sorting the BUS file:

```
echinus count --bus output.correct.sort.bus --ecmap matrix.ec --txnames transcripts.txt --genemap transcripts_to_genes.txt --output echinus_dir --cells barcodes_whitelist.txt
```

The output of echinus is represented in 3 files saved to the directory provided:

1. `matrix.mtx` containing a sparse matrix of gene counts (barcodes in rows),

2. `barcodes.tsv` containing a list of barcodes,

3. `genes.tsv` with tab-separated gene IDs and gene names.

The matrix can be easily loaded into R with [readMM](https://stat.ethz.ch/R-manual/R-devel/library/Matrix/html/externalFormats.html), or with [Read10X](https://www.rdocumentation.org/packages/Seurat/versions/3.1.1/topics/Read10X) when using Seurat for downstream analysis, and into Python with [scipy.io.mmread](https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.mmread.html) or [scanpy.api.read_10x_mtx](https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.api.read_10x_mtx.html).

# Implementation details

Echinus uses [cobra](https://github.com/spf13/cobra)'s implementation of CLI.


