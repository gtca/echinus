![Echinus esculentus](img/cover.png)

# echinus

Get a gene expression table from counts for equivalence classes (ECs -> GE table)

Echinus provides a fast & easy way to summarise transcript compatibility counts to gene abundances. Counts for equivalence classes corresponding to one gene are summed up in order to create a cells x genes matrix commonly used for downstream analysis.

# Usage

```go
echinus -tcc matrix.tcc.mtx -ecmap matrix.ec -txnames transcripts.txt -genemap transcript_to_gene.txt -output echinus_dir/
```
