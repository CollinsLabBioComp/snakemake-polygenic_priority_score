# Parameters

- `study`: Dictionary of GWAS summary statistics with `key`: `<file>`. PoPS will be performed for every GWAS
- `features`: list of feature files to include. The first column of each file should be ensembl gene IDs with header `ENSGID`
- `gene_locations`: file describing genomic locations of genes with the following columns: `ENSGID`, `NAME`, `CHR`, `START`, `END`, `TSS`
- `chromosomes`: list of chromosomes to retain
- `loco`: \[`true`|`false`\]. Whether to employ leave-one-chromosome-out strategy during PoPS analysis (see [PoPS](https://github.com/FinucaneLab/pops))

## Stage 1: MAGMA

- `annotate`: gene annotation setup
  - `run_process`: if false, uses annotation file expected at `data/magma.genes.annot`
  - `window_upstream_kb`: Uses variants within upstream window (in kb)
  - `window_downstream_kb`: Uses variants within downstream window (in kb)
- `score_genes`: gene scoring setup
  - `ld_reference`: BED file for LD reference
  - `sumstat_n_col`: Column in `study` files with sample size (**assumes** same format for all GWAS)
  - `sumstat_variant_id_col`: Column in `study` files with variant IDs (**assumes** same format for all GWAS)
  - `sumstat_pvalue_col`: Column in `study` files with p-values (**assumes** same format for all GWAS)

## Stage 2: PoPS

- `control_features`: file containing control features. Every feauture should be on its own line. Every feature should be found in one of the files described in `features`
- `permute_controls`: \[`true`|`false`\]. Permute controls as well when calculating empirical p-values
- `keep_hla_genes`: \[`true`|`false`\]. Keep HLA genes

## Stage 3: Empricial p-values

- `n_iterations`: number of null permutations
- `grouped`: \[`true`|`false`\]. If `true`, calculates p-values on a per-gene basis. If `false`, uses all genes to calculate p-values.
- `use_absolute_value`: \[`true`|`false`\]. Use absolute value of PoPS scores when calculating p-values. 
- `batch_iterations`: \[`true`|`false`\]. Batch iterations into single jobs to reduce burden on HPCs. 
- `n_batch`: number of batches to use when submitting jobs. Every batch will run `n_iterations` / `n_batch` iterations