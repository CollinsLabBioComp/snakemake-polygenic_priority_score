
# Output directory structure

Below is the structure of the results directory. Pathways will be populated with parameters set in the [params file](../config_analysis.yml). An example of a parameters file is found in the [demo params](../demo/config_analysis.yml).

```bash
results
├── \[files: shared files (i.e., MAGMA, feature prep files)\]
├── gwas_1
│   ├── empirical_pops_pvals.tsv: file containing PoPS scores and empirical p-values
│   ├── empirical_pops_pvals__annotated.tsv: file containing PoPS scores, empirical p-values, and gene information
│   ├── \[files: MAGMA intermediate files\]
│   ├── plots
│   ... etc. ...
├── gwas_2
```
