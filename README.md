# parse_cuffdiff
Report HIDATA, expressed, and significant genes from the output of Cuffdiff.

[Cuffdiff](http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/), which is part of the [Cufflinks toolkit](http://cole-trapnell-lab.github.io/cufflinks/), is used in Bioinformatics to compare expression levels of genes and transcripts in RNA-Seq experiments.

parse_cuffdiff.py takes the output of Cuffdiff and reports several things:

* the HIDATA genes
* the expressed genes
* the significant genes

## Input

Cuffdiff produces several files, three of which are required to run parse_cuffdiff.py:

* gene_exp.diff
* genes.fpkm_tracking
* genes.read_group_tracking

## Output

Four tab-separated value (TSV) text files are produced by parse_cuffdiff.py:

* hidata_genes.tsv - HIDATA genes
* expressed_genes.tsv - expressed genes
* significant_genes.tsv - signficant genes
* gct.tsv - a GCT file

More information about the GCT file format can be found at [The Broad Institute](https://software.broadinstitute.org/software/igv/GCT).

## Requirements

1. [Python3](https://www.python.org/) - Programming Language
2. [Pandas](https://pandas.pydata.org/) - Data Analysis Library

# Usage

python3 ./parse_cuffdiff.py

## Terminal Output

The following example illustrates the output written to the terminal while parse_cuffdiff.py is running:

```
***
* parse_cuffdiff.py
* Parse Cuffdiff output into several files.
***

Input file: genes.fpkm_tracking
Rows: 58024
Columns: 33
Empty cells: 0

Input file: genes.read_group_tracking
Rows: 1044432
Columns: 9
Empty cells: 0

Input file: gene_exp.diff
Rows: 870360
Columns: 14
Empty cells: 566

gene_exp.diff significant=yes...
Rows: 54403
Columns: 14
Empty cells: 221

HIDATA genes...
Columns after removal: 27
Rows with HIDATA: 12
Output file: hidata_genes.tsv

Pivot table & status=OK & mean column...
Rows: 58023
Columns: 19
Empty cells: 132

Expressed genes (mean > 0)...
Rows: 28716
Columns: 19
Empty cells: 65
Output file: expressed_genes.tsv

Significant genes...
Rows: 10336
Columns: 20
Empty cells: 23
Output file: significant_genes.tsv

Making GCT file...
Rows: 10339
Columns: 20
Empty cells: 23
Output file: gct.tsv

Done!
```
