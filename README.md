# XCIE-EM

X-chromosome inactivation escape and phasing from single cell sequencing data using an expectation-maximization algorithm.

### Compiling

```
git clone https://github.com/maickrau/XCIE-EM.git
cd XCIE-EM
git submodule update --init --recursive
make all
```

### Running

First run [scReadCounts](https://horvathlab.github.io/NGS/SCReadCounts/) to count the ref/alt coverage per variant per cell. Then run:

```
XCIE-EM --input-screadcounts screadcounts.output.tsv --output-prefix EM_output --exclude-PAR --exclude-XIST-grch38
```

This will estimate the variant phases and active chrX per cell. The output is written to files which start with `EM_output`. See the section [output](#output) for explanations of all the output files.

The parameters `--exclude-PAR` `--exclude-XIST-grch38` are not mandatory but they improve the results by removing the PAR and XIST regions, which are difficult to phase from chrX inactivation. For other references you can also use `--exclude-XIST-grch37` or `--exclude-XIST-chm13` or manually with `--exclude-region`.

### Output

- `output.variants.tsv`: Description of variant phases, and phasing confidence.
- `output.cells.tsv`: Description of cell activity (which chrX is active per cell), and activity confidence. The cell activity confidence is estimated only from variants which are estimated to be non-escape.
- `output.cells.withescapevariants.tsv`: Description of cell activity (which chrX is active per cell), and activity confidence. The cell activity confidence is estimated from all variants, which biases the cell activity confidence estimate downwards for cells which have high expression in escape variants. Not recommended to be used unless you know what you are doing.
- `output.pseudobulk.variants.confidence2.tsv`: The pseudobulk expression of each variant. The expression is split by chromosome, and further by whether the chromosome is the active chrX or inactive chrX. Only variants and cells which are confidently phased are included.
- `output.pseudobulk.variants.confidence2.tsv`: The pseudobulk expression of each variant. The expression is split by chromosome, and further by whether the chromosome is the active chrX or inactive chrX. All variants and cells are included, even those whose phasing is uncertain.
- `output.preprocessed_matches.tsv`: A file with the cell vs variant matches preprocessed into a suitable format.

### Parameters

- `--input-screadcounts`: Take the counts from scReadCounts as input. Recommended.
- `--input-preprocessed-table`: Take the counts from a preprocessed table.
- `-o`: Output prefix.
- `--force-phase`: Use variant phasing from trio data to force the phase for some variants. This will also help guide the phasing of other variants.
- `--EM-noise-magnitude`: Noise added to the EM algorithm to avoid local minima. Default 20 usually works well.
- `--EM-noise-decay`: Decay speed of noise in EM algorithm. Default 0.95 usually works well.
- `--EM-random-seed`: Random seed used for EM initialization.
- `--EM-num-runs`: Number of initializations used for EM. The best result is picked out of all runs.
- `--exclude-regions`: Exclude some regions. The reason to exclude regions is that some regions (eg XIST) have higher expression from the inactive chrX, which violates the fundamental assumptions used in the EM algorithm.
- `--exclude-PAR`: Exclude the PAR region. Recommended since the PAR region usually has low phasing accuracy.
- `--exclude-XIST-grch38 --exclude-XIST-grch37 --exclude-XIST-chm13`: Exclude the XIST and TSIX genes. Recommended since XIST has higher expression from the inactive chrX. You should pick the correct one based on your reference.

